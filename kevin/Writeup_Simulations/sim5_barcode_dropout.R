# ==============================================================================
# sim5_barcode_dropout.R
#
# GOAL: Quantify the effect of incomplete barcode overlap and cell subsampling
# on (1) CYFER fate potential estimation, and (2) the Gini coefficient of
# the estimated fate potential distribution.
#
# REVIEWER CONCERNS:
#   Comment 4  — "how stable barcode overlap is across timepoints and how clone
#                size affects CYFER's estimation."
#   Comment 5  — "In the cancer dataset, what is the fraction of barcodes that
#                overlap between baseline, day 10, and week 5?"
#   Comment 6  — "How does sampling bias (especially in smaller clones) affect
#                the estimation of the Fate Potential or Gini Coefficient?
#                Perhaps a simulation experiment akin to the one presented in
#                Figure 3 would be helpful here."
#
# DESIGN:
#   Starting from a complete simulated experiment, we introduce two types
#   of incomplete sampling and measure the degradation in CYFER estimates:
#
#   PART A — Clone-level dropout (t2 barcode overlap):
#     A fraction (1 - p_overlap) of clones are undetected at t2.
#     These clones have 0 observed count in lineage_future_count.
#     Two strategies tested:
#       (i)  Random dropout: any clone may be undetected.
#       (ii) Size-biased dropout: small clones at t2 are more likely to be
#            missed (reflecting lower sequencing depth for rare barcodes).
#     Metric: correlation of CYFER-estimated fate potential vs. true, and
#             Gini coefficient of estimated fate potential.
#
#   PART B — Cell-level subsampling (incomplete cell capture at t1):
#     At t1 (the early time point), only a fraction (p_cell_capture) of cells
#     are sequenced. CYFER is fit using only the subsampled cells.
#     Metric: correlation, Gini coefficient.
#
#   PART C — Clone size effect on Gini estimation:
#     Fix p_overlap = 1.0, vary minimum clone size threshold used when
#     computing the Gini coefficient of fate potentials.
#     Demonstrates how the minimum-size cutoff affects the Gini estimate.
#
# EXPECTED RESULTS:
#   CYFER is robust to moderate dropout because it aggregates over all cells
#   within a clone, and the loss of a few clones from lineage_future_count
#   only slightly reduces the training signal. The Gini coefficient, computed
#   on imputed fate potentials, is more sensitive to cell-level subsampling.
# ==============================================================================

rm(list=ls())
library(multiomeFate)
library(MASS)

set.seed(42)

# ---- Helper functions --------------------------------------------------------

gini_coef <- function(x) {
  x <- pmax(x, 0); x <- sort(x); n <- length(x)
  if (n == 0 || sum(x) == 0) return(0)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

fit_cyfer <- function(X, clone_labels, lineage_future_count) {
  valid_clones <- names(lineage_future_count[lineage_future_count > 0])
  if (length(valid_clones) < 3) return(NULL)
  keep_idx  <- which(clone_labels %in% valid_clones)
  X_sub     <- X[keep_idx, , drop = FALSE]
  clone_sub <- clone_labels[keep_idx]
  lfc_sub   <- lineage_future_count[valid_clones]
  clone_sizes <- table(clone_sub)
  valid_clones2 <- names(clone_sizes[clone_sizes >= 2])
  if (length(valid_clones2) < 3) return(NULL)
  keep2 <- which(clone_sub %in% valid_clones2)
  X_sub <- X_sub[keep2, , drop = FALSE]; clone_sub <- clone_sub[keep2]
  lfc_sub <- lfc_sub[valid_clones2]
  n_folds <- min(3, length(unique(clone_sub)) - 1)
  if (n_folds < 2) return(NULL)
  tryCatch({
    fit_res <- cyfer(cell_features = X_sub, cell_lineage = clone_sub,
                     lineage_future_count = lfc_sub, lambda_initial = 1,
                     lambda_sequence_length = 10, num_folds = n_folds, verbose = 0)
    cyfer_finalize(cell_features = X_sub, cell_lineage = clone_sub,
                   fit_res = fit_res, lineage_future_count = lfc_sub)
  }, error = function(e) NULL)
}

# ---- Baseline simulation parameters ------------------------------------------

n_clones     <- 100
n_per_clone  <- 25
n_cells      <- n_clones * n_per_clone
n_features   <- 20
n_causal     <- 4
sigma_between <- 1.0
sigma_within  <- 0.4
n_replicates  <- 10

true_beta      <- c(rep(1.5, n_causal), rep(0, n_features - n_causal))
feat_names     <- paste0("f", seq_len(n_features))
names(true_beta) <- feat_names

clone_ids    <- rep(seq_len(n_clones), each = n_per_clone)
clone_labels <- paste0("clone:", clone_ids)

# ---- Generate one complete baseline dataset ----------------------------------

generate_base_data <- function(seed_val) {
  set.seed(seed_val)
  cc <- mvrnorm(n_clones, rep(0, n_features), sigma_between^2 * diag(n_features))
  colnames(cc) <- feat_names
  X <- t(sapply(seq_len(n_cells), function(i) {
    cc[clone_ids[i], ] + mvrnorm(1, rep(0, n_features), sigma_within^2 * diag(n_features))
  }))
  rownames(X) <- paste0("cell:", seq_len(n_cells)); colnames(X) <- feat_names
  X <- scale(X)
  true_Z <- as.numeric(X %*% true_beta); names(true_Z) <- rownames(X)
  true_Z <- true_Z - mean(true_Z)
  cell_future <- rpois(n_cells, exp(true_Z)); names(cell_future) <- rownames(X)
  lfc_full <- tapply(cell_future, clone_labels, sum)
  list(X = X, true_Z = true_Z, lfc_full = lfc_full)
}

# Gini of the true fate potential
true_gini <- function(true_Z) gini_coef(exp(true_Z))

# ---- PART A: Clone-level dropout (barcode overlap) ---------------------------

p_overlap_values <- c(0.10, 0.25, 0.50, 0.75, 1.00)

message("=== PART A: Clone-level dropout ===")

results_A_random <- lapply(p_overlap_values, function(p_ov) {
  message("  p_overlap = ", p_ov)
  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    print(paste0("Replicate: ", rep_idx))
    base <- generate_base_data(rep_idx * 11 + round(p_ov * 100))
    X       <- base$X
    true_Z  <- base$true_Z
    lfc_full <- base$lfc_full

    # Random dropout: keep only p_ov fraction of clones
    n_keep   <- max(3, round(p_ov * n_clones))
    kept_clones  <- sample(names(lfc_full), n_keep)
    lfc_dropout  <- lfc_full[kept_clones]
    lfc_dropout  <- lfc_dropout[lfc_dropout > 0]

    final_fit <- fit_cyfer(X, clone_labels, lfc_dropout)
    if (is.null(final_fit)) return(NULL)

    Z_hat <- as.numeric(X %*% final_fit$coefficient_vec[-1]) + final_fit$coefficient_vec[1]
    names(Z_hat) <- rownames(X)

    list(
      cor_z      = cor(true_Z, Z_hat),
      gini_true  = gini_coef(exp(true_Z)),
      gini_hat   = gini_coef(exp(Z_hat)),
      gini_error = abs(gini_coef(exp(Z_hat)) - gini_coef(exp(true_Z)))
    )
  })
  reps <- Filter(Negate(is.null), reps)
  data.frame(
    p_overlap         = p_ov,
    dropout_type      = "random",
    mean_cor_z        = round(mean(sapply(reps, `[[`, "cor_z"),      na.rm = TRUE), 3),
    sd_cor_z          = round(sd(  sapply(reps, `[[`, "cor_z"),      na.rm = TRUE), 3),
    mean_gini_true    = round(mean(sapply(reps, `[[`, "gini_true"),  na.rm = TRUE), 3),
    mean_gini_hat     = round(mean(sapply(reps, `[[`, "gini_hat"),   na.rm = TRUE), 3),
    mean_gini_error   = round(mean(sapply(reps, `[[`, "gini_error"), na.rm = TRUE), 3),
    n_converged       = length(reps),
    stringsAsFactors  = FALSE
  )
})

# Size-biased dropout: small clones (few cells at t2) more likely to be missed
results_A_biased <- lapply(p_overlap_values, function(p_ov) {
  message("  p_overlap = ", p_ov, " (size-biased)")
  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    print(paste0("Replicate: ", rep_idx))
    set.seed(rep_idx * 11 + round(p_ov * 100) + 5000)
    base <- generate_base_data(rep_idx * 11 + round(p_ov * 100))
    lfc_full <- base$lfc_full; X <- base$X; true_Z <- base$true_Z

    # Size-biased: probability of detection proportional to clone size
    detection_prob <- pmin(lfc_full / quantile(lfc_full, 0.5), 1) * p_ov
    detected <- rbinom(length(lfc_full), 1, prob = detection_prob) == 1
    lfc_dropout <- lfc_full[detected & lfc_full > 0]
    if (length(lfc_dropout) < 3) return(NULL)

    final_fit <- fit_cyfer(X, clone_labels, lfc_dropout)
    if (is.null(final_fit)) return(NULL)

    Z_hat <- as.numeric(X %*% final_fit$coefficient_vec[-1]) + final_fit$coefficient_vec[1]
    names(Z_hat) <- rownames(X)

    list(cor_z = cor(true_Z, Z_hat),
         gini_true  = gini_coef(exp(true_Z)),
         gini_hat   = gini_coef(exp(Z_hat)),
         gini_error = abs(gini_coef(exp(Z_hat)) - gini_coef(exp(true_Z))))
  })
  reps <- Filter(Negate(is.null), reps)
  data.frame(
    p_overlap       = p_ov, dropout_type = "size_biased",
    mean_cor_z      = round(mean(sapply(reps, `[[`, "cor_z"),      na.rm = TRUE), 3),
    sd_cor_z        = round(sd(  sapply(reps, `[[`, "cor_z"),      na.rm = TRUE), 3),
    mean_gini_true  = round(mean(sapply(reps, `[[`, "gini_true"),  na.rm = TRUE), 3),
    mean_gini_hat   = round(mean(sapply(reps, `[[`, "gini_hat"),   na.rm = TRUE), 3),
    mean_gini_error = round(mean(sapply(reps, `[[`, "gini_error"), na.rm = TRUE), 3),
    n_converged     = length(reps), stringsAsFactors = FALSE
  )
})

results_A_df <- rbind(do.call(rbind, results_A_random),
                      do.call(rbind, results_A_biased))
cat("\n--- Part A results (clone dropout) ---\n")
print(results_A_df)

# ---- PART B: Cell-level subsampling at t1 ------------------------------------

p_cell_values <- c(0.10, 0.25, 0.50, 0.75, 1.00)

message("\n=== PART B: Cell-level subsampling at t1 ===")

results_B <- lapply(p_cell_values, function(p_cell) {
  message("  p_cell_capture = ", p_cell)
  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    print(paste0("Replicate: ", rep_idx))
    base <- generate_base_data(rep_idx * 7 + round(p_cell * 100))
    X <- base$X; true_Z <- base$true_Z; lfc_full <- base$lfc_full

    # Subsample cells (at the t1 collection step)
    n_obs      <- max(10, round(p_cell * n_cells))
    obs_idx    <- sample(n_cells, n_obs, replace = FALSE)
    X_obs      <- X[obs_idx, , drop = FALSE]
    clone_obs  <- clone_labels[obs_idx]
    Z_obs      <- true_Z[obs_idx]

    # Only keep clones with >= 2 observed cells
    clone_counts <- table(clone_obs)
    valid_cl     <- names(clone_counts[clone_counts >= 2])
    lfc_obs      <- lfc_full[valid_cl[valid_cl %in% names(lfc_full)]]
    lfc_obs      <- lfc_obs[lfc_obs > 0]
    if (length(lfc_obs) < 3) return(NULL)

    final_fit <- fit_cyfer(X_obs, clone_obs, lfc_obs)
    if (is.null(final_fit)) return(NULL)

    Z_hat <- as.numeric(X_obs %*% final_fit$coefficient_vec[-1]) +
             final_fit$coefficient_vec[1]
    names(Z_hat) <- rownames(X_obs)

    list(cor_z      = cor(Z_obs, Z_hat),
         gini_true  = gini_coef(exp(Z_obs)),
         gini_hat   = gini_coef(exp(Z_hat)),
         gini_error = abs(gini_coef(exp(Z_hat)) - gini_coef(exp(Z_obs))))
  })
  reps <- Filter(Negate(is.null), reps)
  data.frame(
    p_cell_capture  = p_cell,
    mean_cor_z      = round(mean(sapply(reps, `[[`, "cor_z"),      na.rm = TRUE), 3),
    sd_cor_z        = round(sd(  sapply(reps, `[[`, "cor_z"),      na.rm = TRUE), 3),
    mean_gini_true  = round(mean(sapply(reps, `[[`, "gini_true"),  na.rm = TRUE), 3),
    mean_gini_hat   = round(mean(sapply(reps, `[[`, "gini_hat"),   na.rm = TRUE), 3),
    mean_gini_error = round(mean(sapply(reps, `[[`, "gini_error"), na.rm = TRUE), 3),
    n_converged     = length(reps), stringsAsFactors = FALSE
  )
})
results_B_df <- do.call(rbind, results_B)

cat("\n--- Part B results (cell subsampling at t1) ---\n")
print(results_B_df)

# ---- PART C: Minimum clone size threshold for Gini ---------------------------

message("\n=== PART C: Gini coefficient vs minimum clone size threshold ===")

min_clone_sizes <- c(1, 5, 10, 15, 20, 30)

set.seed(777)
base_main <- generate_base_data(777)
X_main    <- base_main$X; true_Z_main <- base_main$true_Z; lfc_main <- base_main$lfc_full

final_fit_main <- fit_cyfer(X_main, clone_labels, lfc_main)
if (!is.null(final_fit_main)) {
  Z_hat_main <- as.numeric(X_main %*% final_fit_main$coefficient_vec[-1]) +
               final_fit_main$coefficient_vec[1]
  names(Z_hat_main) <- rownames(X_main)

  gini_vs_threshold <- lapply(min_clone_sizes, function(min_sz) {
    # Restrict to clones with at least min_sz cells
    clone_size_vec <- table(clone_labels)
    large_clones   <- names(clone_size_vec[clone_size_vec >= min_sz])
    keep           <- which(clone_labels %in% large_clones)

    gini_true <- gini_coef(exp(true_Z_main[keep]))
    gini_hat  <- gini_coef(exp(Z_hat_main[keep]))

    data.frame(min_clone_size = min_sz,
               n_clones_used  = length(large_clones),
               n_cells_used   = length(keep),
               gini_true      = round(gini_true, 4),
               gini_hat       = round(gini_hat,  4),
               gini_error     = round(abs(gini_hat - gini_true), 4),
               stringsAsFactors = FALSE)
  })
  results_C_df <- do.call(rbind, gini_vs_threshold)
  cat("\n--- Part C results (Gini vs clone size threshold) ---\n")
  print(results_C_df)
} else {
  results_C_df <- NULL
  cat("Part C: could not fit CYFER on main dataset\n")
}

# ---- Save results ------------------------------------------------------------

filepath <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup_Simulations/"
saveRDS(
  list(partA = results_A_df,
       partB = results_B_df,
       partC = results_C_df,
       params = list(n_clones = n_clones, n_per_clone = n_per_clone,
                     n_features = n_features, n_causal = n_causal,
                     n_replicates = n_replicates)),
  file = paste0(filepath, "sim5_barcode_dropout_results.rds")
)

message("sim5 complete. Results saved to sim5_barcode_dropout_results.rds")
