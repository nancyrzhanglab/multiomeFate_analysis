# ==============================================================================
# sim4_power_analysis.R
#
# GOAL: Determine the minimum experimental requirements for reliable CYFER
# estimation, as a function of the number of clones, cells per clone, and
# barcode capture efficacy.
#
# REVIEWER CONCERNS:
#   Comment 1d — "(d) different cell numbers and barcoding efficacies (to give
#                some estimates of power)."
#   Comment 4  — "how stable barcode overlap is across timepoints and how clone
#                size affects CYFER's estimation."
#   Comment 5  — "barcode statistics."
#   Comment 9  — "further simulation to determine the minimum number of clones
#                for stable performance."
#
# DESIGN:
#   Power is assessed along three independent axes:
#
#   AXIS 1 — Number of clones (K):
#     Fix cells per clone = 20, vary K in {10, 25, 50, 100, 200, 500}.
#     Measures: How many clones are needed for stable coefficient estimation?
#
#   AXIS 2 — Cells per clone (n_per_clone):
#     Fix K = 75, vary n_per_clone in {3, 5, 10, 20, 50}.
#     Measures: How many cells per clone are needed?
#
#   AXIS 3 — Barcode capture efficacy (p_capture):
#     Fix K = 75, n_per_clone = 20.
#     At t2, only a fraction p_capture of clones are detected.
#     Simulate by randomly removing clones from lineage_future_count.
#     Vary p_capture in {0.10, 0.25, 0.50, 0.75, 1.00}.
#     Measures: How much barcode dropout can CYFER tolerate?
#
#   For each parameter combination, n_replicates=25 independent datasets are
#   generated and CYFER is fit to each. Performance is measured as:
#     - Pearson correlation of estimated vs. true cell fate potential.
#     - Jaccard index of identified feature set.
#     - Coefficient stability (mean absolute deviation of estimated beta from
#       true beta, normalized by true signal magnitude).
#
# EXPECTED RESULTS:
#   Performance should increase monotonically with K, n_per_clone, and
#   p_capture, but with diminishing returns. The simulation should identify
#   practical thresholds below which CYFER becomes unreliable.
# ==============================================================================

library(multiomeFate)
library(MASS)

set.seed(42)

# ---- Helper functions --------------------------------------------------------

gini_coef <- function(x) {
  x <- pmax(x, 0); x <- sort(x); n <- length(x)
  if (n == 0 || sum(x) == 0) return(0)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

jaccard <- function(a, b) {
  u <- union(a, b); if (length(u) == 0) return(0)
  length(intersect(a, b)) / length(u)
}

top_features_by_cor <- function(feature_mat, target, k) {
  cors <- apply(feature_mat, 2, function(col) {
    if (sd(col) < 1e-10) return(0); cor(col, target)
  })
  names(sort(abs(cors), decreasing = TRUE))[seq_len(min(k, ncol(feature_mat)))]
}

# ---- Default simulation parameters -------------------------------------------

n_features   <- 20
n_causal     <- 4
sigma_between <- 1.0
sigma_within  <- 0.5

true_beta      <- c(rep(1.2, n_causal), rep(0, n_features - n_causal))
feat_names     <- paste0("f", seq_len(n_features))
names(true_beta) <- feat_names
causal_features <- feat_names[seq_len(n_causal)]

n_replicates <- 25

# ---- Core CYFER fitting function ---------------------------------------------

fit_cyfer_and_eval <- function(X, clone_labels, lineage_future_count, true_Z) {
  valid_clones <- names(lineage_future_count[lineage_future_count > 0])
  if (length(valid_clones) < 3) return(NULL)

  keep_idx  <- which(clone_labels %in% valid_clones)
  X_sub     <- X[keep_idx, , drop = FALSE]
  Z_sub     <- true_Z[keep_idx]
  clone_sub <- clone_labels[keep_idx]
  lfc_sub   <- lineage_future_count[valid_clones]

  # Need at least 3 cells in each clone for stable CV
  clone_sizes <- table(clone_sub)
  valid_clones2 <- names(clone_sizes[clone_sizes >= 2])
  if (length(valid_clones2) < 3) return(NULL)

  keep_idx2 <- which(clone_sub %in% valid_clones2)
  X_sub     <- X_sub[keep_idx2, , drop = FALSE]
  Z_sub     <- Z_sub[keep_idx2]
  clone_sub <- clone_sub[keep_idx2]
  lfc_sub   <- lfc_sub[valid_clones2]

  n_folds <- min(3, length(unique(clone_sub)) - 1)
  if (n_folds < 2) return(NULL)

  tryCatch({
    fit_res <- cyfer(
      cell_features        = X_sub,
      cell_lineage         = clone_sub,
      lineage_future_count = lfc_sub,
      lambda_initial       = 1,
      lambda_sequence_length = 10,
      num_folds            = n_folds,
      verbose              = 0
    )
    final_fit <- cyfer_finalize(
      cell_features        = X_sub,
      cell_lineage         = clone_sub,
      fit_res              = fit_res,
      lineage_future_count = lfc_sub
    )
    Z_hat     <- as.numeric(X_sub %*% final_fit$coefficient_vec[-1]) +
                 final_fit$coefficient_vec[1]
    names(Z_hat) <- rownames(X_sub)

    cor_z     <- cor(Z_sub, Z_hat)
    est_feat  <- top_features_by_cor(X_sub, Z_hat, k = n_causal)
    jac_val   <- jaccard(causal_features, est_feat)

    # Coefficient recovery (alignment between true_beta and estimated beta)
    coef_est  <- final_fit$coefficient_vec[-1]  # drop intercept
    if (all(feat_names %in% names(coef_est))) {
      coef_est <- coef_est[feat_names]
    }
    coef_cor  <- cor(true_beta, coef_est)

    list(cor_z = cor_z, jaccard = jac_val, coef_cor = coef_cor,
         n_valid_clones = length(valid_clones2),
         n_valid_cells  = nrow(X_sub))
  }, error = function(e) NULL)
}

# ---- AXIS 1: Number of clones ------------------------------------------------

K_values     <- c(10, 25, 50, 100, 200, 500)
n_per_clone  <- 20   # fixed cells per clone

message("=== AXIS 1: Varying number of clones (n_per_clone = ", n_per_clone, ") ===")

results_K <- lapply(K_values, function(K) {
  message("  K = ", K)
  n_cells  <- K * n_per_clone
  clone_ids    <- rep(seq_len(K), each = n_per_clone)
  clone_labels <- paste0("clone:", clone_ids)

  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    set.seed(rep_idx * 13 + K)

    cc <- mvrnorm(K, mu = rep(0, n_features), Sigma = sigma_between^2 * diag(n_features))
    colnames(cc) <- feat_names

    X <- t(sapply(seq_len(n_cells), function(i) {
      cc[clone_ids[i], ] + mvrnorm(1, rep(0, n_features), sigma_within^2 * diag(n_features))
    }))
    rownames(X) <- paste0("cell:", seq_len(n_cells)); colnames(X) <- feat_names
    X <- scale(X)

    true_Z <- as.numeric(X %*% true_beta); names(true_Z) <- rownames(X)
    true_Z <- true_Z - mean(true_Z)

    cell_future <- rpois(n_cells, exp(true_Z))
    names(cell_future) <- rownames(X)
    lineage_future_count <- tapply(cell_future, clone_labels, sum)

    fit_cyfer_and_eval(X, clone_labels, lineage_future_count, true_Z)
  })
  reps <- Filter(Negate(is.null), reps)

  data.frame(
    K                  = K,
    mean_cor_z         = round(mean(sapply(reps, `[[`, "cor_z"),    na.rm = TRUE), 3),
    sd_cor_z           = round(sd(  sapply(reps, `[[`, "cor_z"),    na.rm = TRUE), 3),
    mean_jaccard       = round(mean(sapply(reps, `[[`, "jaccard"),  na.rm = TRUE), 3),
    mean_coef_cor      = round(mean(sapply(reps, `[[`, "coef_cor"), na.rm = TRUE), 3),
    n_converged        = length(reps),
    stringsAsFactors   = FALSE
  )
})
results_K_df <- do.call(rbind, results_K)

cat("\n--- Axis 1 results (varying number of clones) ---\n")
print(results_K_df)

# ---- AXIS 2: Cells per clone -------------------------------------------------

n_per_clone_values <- c(3, 5, 10, 20, 50)
K_fixed            <- 75

message("\n=== AXIS 2: Varying cells per clone (K = ", K_fixed, ") ===")

results_npc <- lapply(n_per_clone_values, function(npc) {
  message("  n_per_clone = ", npc)
  n_cells      <- K_fixed * npc
  clone_ids    <- rep(seq_len(K_fixed), each = npc)
  clone_labels <- paste0("clone:", clone_ids)

  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    set.seed(rep_idx * 19 + npc * 3)

    cc <- mvrnorm(K_fixed, mu = rep(0, n_features), Sigma = sigma_between^2 * diag(n_features))
    colnames(cc) <- feat_names

    X <- t(sapply(seq_len(n_cells), function(i) {
      cc[clone_ids[i], ] + mvrnorm(1, rep(0, n_features), sigma_within^2 * diag(n_features))
    }))
    rownames(X) <- paste0("cell:", seq_len(n_cells)); colnames(X) <- feat_names
    X <- scale(X)

    true_Z <- as.numeric(X %*% true_beta); names(true_Z) <- rownames(X)
    true_Z <- true_Z - mean(true_Z)

    cell_future <- rpois(n_cells, exp(true_Z)); names(cell_future) <- rownames(X)
    lineage_future_count <- tapply(cell_future, clone_labels, sum)

    fit_cyfer_and_eval(X, clone_labels, lineage_future_count, true_Z)
  })
  reps <- Filter(Negate(is.null), reps)

  data.frame(
    n_per_clone        = npc,
    mean_cor_z         = round(mean(sapply(reps, `[[`, "cor_z"),    na.rm = TRUE), 3),
    sd_cor_z           = round(sd(  sapply(reps, `[[`, "cor_z"),    na.rm = TRUE), 3),
    mean_jaccard       = round(mean(sapply(reps, `[[`, "jaccard"),  na.rm = TRUE), 3),
    mean_coef_cor      = round(mean(sapply(reps, `[[`, "coef_cor"), na.rm = TRUE), 3),
    n_converged        = length(reps),
    stringsAsFactors   = FALSE
  )
})
results_npc_df <- do.call(rbind, results_npc)

cat("\n--- Axis 2 results (varying cells per clone) ---\n")
print(results_npc_df)

# ---- AXIS 3: Barcode capture efficacy ----------------------------------------

p_capture_values <- c(0.10, 0.25, 0.50, 0.75, 1.00)
K_axis3          <- 75
npc_axis3        <- 20
n_cells_axis3    <- K_axis3 * npc_axis3
clone_ids_axis3  <- rep(seq_len(K_axis3), each = npc_axis3)
clone_labels_axis3 <- paste0("clone:", clone_ids_axis3)

message("\n=== AXIS 3: Varying barcode capture efficacy (K=", K_axis3,
        ", n_per_clone=", npc_axis3, ") ===")

results_cap <- lapply(p_capture_values, function(p_cap) {
  message("  p_capture = ", p_cap)

  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    set.seed(rep_idx * 23 + round(p_cap * 100))

    cc <- mvrnorm(K_axis3, mu = rep(0, n_features), Sigma = sigma_between^2 * diag(n_features))
    colnames(cc) <- feat_names

    X <- t(sapply(seq_len(n_cells_axis3), function(i) {
      cc[clone_ids_axis3[i], ] +
        mvrnorm(1, rep(0, n_features), sigma_within^2 * diag(n_features))
    }))
    rownames(X) <- paste0("cell:", seq_len(n_cells_axis3)); colnames(X) <- feat_names
    X <- scale(X)

    true_Z <- as.numeric(X %*% true_beta); names(true_Z) <- rownames(X)
    true_Z <- true_Z - mean(true_Z)

    # Full future counts
    cell_future <- rpois(n_cells_axis3, exp(true_Z)); names(cell_future) <- rownames(X)
    lineage_future_count_full <- tapply(cell_future, clone_labels_axis3, sum)

    # Simulate barcode dropout: only p_cap fraction of clones observed at t2
    all_clones     <- names(lineage_future_count_full)
    n_observed     <- max(3, round(p_cap * length(all_clones)))
    observed_clones <- sample(all_clones, n_observed, replace = FALSE)
    lineage_future_count <- lineage_future_count_full[observed_clones]
    lineage_future_count <- lineage_future_count[lineage_future_count > 0]

    if (length(lineage_future_count) < 3) return(NULL)

    fit_cyfer_and_eval(X, clone_labels_axis3, lineage_future_count, true_Z)
  })
  reps <- Filter(Negate(is.null), reps)

  data.frame(
    p_capture          = p_cap,
    mean_cor_z         = round(mean(sapply(reps, `[[`, "cor_z"),    na.rm = TRUE), 3),
    sd_cor_z           = round(sd(  sapply(reps, `[[`, "cor_z"),    na.rm = TRUE), 3),
    mean_jaccard       = round(mean(sapply(reps, `[[`, "jaccard"),  na.rm = TRUE), 3),
    mean_coef_cor      = round(mean(sapply(reps, `[[`, "coef_cor"), na.rm = TRUE), 3),
    n_converged        = length(reps),
    stringsAsFactors   = FALSE
  )
})
results_cap_df <- do.call(rbind, results_cap)

cat("\n--- Axis 3 results (varying barcode capture efficacy) ---\n")
print(results_cap_df)

# ---- Final summary -----------------------------------------------------------

cat("\n=== POWER ANALYSIS SUMMARY ===\n")
cat("Minimum recommended thresholds (mean_cor_z >= 0.6):\n")
cat("  Clones needed:        K >=",
    results_K_df$K[min(which(results_K_df$mean_cor_z >= 0.6))], "\n")
cat("  Cells per clone:      n_per_clone >=",
    results_npc_df$n_per_clone[min(which(results_npc_df$mean_cor_z >= 0.6))], "\n")
cat("  Barcode capture rate: p_capture >=",
    results_cap_df$p_capture[min(which(results_cap_df$mean_cor_z >= 0.6))], "\n")

# ---- Save results ------------------------------------------------------------

saveRDS(
  list(axis1_clones   = results_K_df,
       axis2_cellsize  = results_npc_df,
       axis3_capture   = results_cap_df,
       params = list(n_features = n_features, n_causal = n_causal,
                     sigma_between = sigma_between, sigma_within = sigma_within,
                     n_replicates = n_replicates)),
  file = "sim4_power_analysis_results.rds"
)

message("sim4 complete. Results saved to sim4_power_analysis_results.rds")
