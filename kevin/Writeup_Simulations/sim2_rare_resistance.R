# ==============================================================================
# sim2_rare_resistance.R
#
# GOAL: Determine CYFER's ability to detect rare "resistant" progenitor cells
# whose fate potential is far higher than the population average, across a
# range of prevalence levels.
#
# REVIEWER CONCERN: Comment 1b — "different distributions of phenotype
# proportions (e.g. when resistance is very rare vs. moderately common)."
# Also relevant to Comment 9 regarding detecting rare features predictive of
# resistance.
#
# DESIGN:
#   - n=2000 cells, K=100 clones, p=30 features.
#   - A fraction f of cells are "resistant": they express a specific gene
#     program (features 1-5 at high level) and have high fate potential.
#   - The remaining (1-f) cells are "non-resistant": low expression of
#     resistance features and low/zero fate potential.
#   - f is varied: 0.01, 0.02, 0.05, 0.10, 0.25, 0.50.
#   - Cells are randomly assigned to clones so most clones contain a mixture.
#   - CYFER is fit and cell-level fate potentials are computed.
#   - Performance: AUROC, sensitivity, specificity, Jaccard for feature
#     identification.
#
# BIOLOGICAL INTERPRETATION:
#   In the cancer resistance context, rare pre-existing cells (1-5% of the
#   population) are sufficient to seed long-term resistance. The question is
#   whether CYFER can identify these rare cells from the lineage data.
#
# EXPECTED RESULT:
#   CYFER should have high sensitivity even when f=0.01 (1% resistant), because
#   the exponential growth model amplifies the signal from rare high-potential
#   cells. Performance degrades when f is very small AND clone sizes are small
#   (some clones may have 0 resistant cells, making the signal undetectable).
# ==============================================================================

library(multiomeFate)
library(MASS)

set.seed(42)

# ---- Helper functions --------------------------------------------------------

gini_coef <- function(x) {
  x <- pmax(x, 0)
  x <- sort(x)
  n <- length(x)
  if (n == 0 || sum(x) == 0) return(0)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(0)
  length(intersect(a, b)) / length(u)
}

# Area under ROC curve (Wilcoxon statistic)
auroc <- function(scores, labels) {
  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)
  if (n_pos == 0 || n_neg == 0) return(0.5)
  wilcox.test(scores[labels == 1], scores[labels == 0],
              alternative = "greater")$statistic / (n_pos * n_neg)
}

top_features_by_cor <- function(feature_mat, target, k) {
  cors <- apply(feature_mat, 2, function(col) {
    if (sd(col) < 1e-10) return(0)
    cor(col, target)
  })
  names(sort(abs(cors), decreasing = TRUE))[seq_len(min(k, ncol(feature_mat)))]
}

# ---- Simulation parameters ---------------------------------------------------

n_cells      <- 2000
n_clones     <- 100
n_features   <- 30
n_causal     <- 5        # features distinguishing resistant from non-resistant

# Resistance fractions to test
f_values     <- c(0.01, 0.02, 0.05, 0.10, 0.25, 0.50)

# Expression parameters
mu_resistant    <- 2.0    # resistant cells have high expression of causal features
mu_background   <- 0.0    # non-resistant cells have low expression
sigma_expr      <- 0.8    # within-group expression noise

# Fate potential parameters
Z_resistant     <-  3.0   # log expected progeny for resistant cells
Z_nonresistant  <- -1.5   # log expected progeny for non-resistant cells
Z_noise         <-  0.3   # cell-level fate noise

n_replicates    <- 25

# Feature names
feat_names     <- paste0("f", seq_len(n_features))
causal_features <- paste0("f", seq_len(n_causal))

# ---- Run simulation for one replicate at one prevalence ----------------------

run_one <- function(f, seed_val) {
  set.seed(seed_val)

  n_resistant    <- max(1, round(f * n_cells))
  n_nonresistant <- n_cells - n_resistant

  # Assign resistance status
  is_resistant <- c(rep(TRUE, n_resistant), rep(FALSE, n_nonresistant))
  cell_names   <- paste0("cell:", seq_len(n_cells))

  # Generate features
  # Causal features: resistant cells have high expression
  causal_mat <- rbind(
    matrix(rnorm(n_resistant    * n_causal, mu_resistant,  sigma_expr),
           nrow = n_resistant),
    matrix(rnorm(n_nonresistant * n_causal, mu_background, sigma_expr),
           nrow = n_nonresistant)
  )

  # Non-causal features: random noise for all cells
  noise_mat <- matrix(rnorm(n_cells * (n_features - n_causal), 0, 1),
                      nrow = n_cells)

  X <- cbind(causal_mat, noise_mat)
  rownames(X) <- cell_names
  colnames(X) <- feat_names
  X <- scale(X)

  # True fate potential
  true_Z <- ifelse(is_resistant, Z_resistant, Z_nonresistant) +
            rnorm(n_cells, 0, Z_noise)
  names(true_Z) <- cell_names

  # Random clone assignment (balanced)
  cells_per_clone <- n_cells %/% n_clones
  clone_ids       <- rep(seq_len(n_clones), each = cells_per_clone)[seq_len(n_cells)]
  # Shuffle assignment so resistant cells spread across clones
  shuffle_idx     <- sample(n_cells)
  clone_ids       <- clone_ids[order(shuffle_idx)]
  clone_labels    <- paste0("clone:", clone_ids)

  # Generate future clone sizes (exponential model)
  cell_future <- rpois(n_cells, lambda = exp(true_Z))
  names(cell_future) <- cell_names
  lineage_future_count <- tapply(cell_future, clone_labels, sum)

  # Separate clones with and without resistant cells
  has_resistant <- tapply(is_resistant, clone_labels, any)

  valid_clones <- names(lineage_future_count[lineage_future_count > 0])
  if (length(valid_clones) < 5) return(NULL)

  keep_idx  <- which(clone_labels %in% valid_clones)
  X_sub     <- X[keep_idx, , drop = FALSE]
  Z_sub     <- true_Z[keep_idx]
  ir_sub    <- is_resistant[keep_idx]
  clone_sub <- clone_labels[keep_idx]
  lfc_sub   <- lineage_future_count[valid_clones]

  tryCatch({
    fit_res <- cyfer(
      cell_features        = X_sub,
      cell_lineage         = clone_sub,
      lineage_future_count = lfc_sub,
      lambda_initial       = 1,
      lambda_sequence_length = 10,
      num_folds            = 3,
      verbose              = 0
    )
    final_fit <- cyfer_finalize(
      cell_features        = X_sub,
      cell_lineage         = clone_sub,
      fit_res              = fit_res,
      lineage_future_count = lfc_sub
    )

    Z_hat <- as.numeric(X_sub %*% final_fit$coefficient_vec[-1]) +
             final_fit$coefficient_vec[1]
    names(Z_hat) <- rownames(X_sub)

    # Cell-level metrics
    auc_val   <- auroc(Z_hat, as.integer(ir_sub))
    cor_z     <- cor(Z_sub, Z_hat)

    # Feature-level metrics (top n_causal*2 features detected)
    est_causal <- top_features_by_cor(X_sub, Z_hat, k = n_causal * 2)
    jac_val    <- jaccard(causal_features,
                          est_causal[seq_len(n_causal)])

    # Sensitivity / specificity at optimal threshold (Youden's J)
    thresholds <- quantile(Z_hat, probs = seq(0.1, 0.9, by = 0.1))
    youden <- sapply(thresholds, function(thr) {
      tp <- sum(Z_hat >= thr &  ir_sub)
      fn <- sum(Z_hat <  thr &  ir_sub)
      tn <- sum(Z_hat <  thr & !ir_sub)
      fp <- sum(Z_hat >= thr & !ir_sub)
      sens <- tp / (tp + fn + 1e-10)
      spec <- tn / (tn + fp + 1e-10)
      sens + spec - 1
    })
    best_thr  <- thresholds[which.max(youden)]
    tp <- sum(Z_hat >= best_thr &  ir_sub)
    fn <- sum(Z_hat <  best_thr &  ir_sub)
    tn <- sum(Z_hat <  best_thr & !ir_sub)
    fp <- sum(Z_hat >= best_thr & !ir_sub)
    sens <- tp / (tp + fn + 1e-10)
    spec <- tn / (tn + fp + 1e-10)

    # Fraction of clones with resistant cells that have high imputed score
    clone_imputed <- tapply(Z_hat, clone_sub, mean)
    clone_has_res <- has_resistant[names(clone_imputed)]
    clone_auc     <- auroc(clone_imputed, as.integer(clone_has_res))

    list(f = f, seed = seed_val, auroc_cell = auc_val,
         auroc_clone = clone_auc, cor_z = cor_z,
         sensitivity = sens, specificity = spec,
         jaccard_features = jac_val,
         n_resistant = n_resistant,
         n_valid_clones = length(valid_clones))
  }, error = function(e) NULL)
}

# ---- Run all combinations ----------------------------------------------------

message("Running ", n_replicates, " replicates x ", length(f_values), " prevalence values...")

all_results <- lapply(f_values, function(f) {
  message("  f = ", f)
  res_list <- lapply(seq_len(n_replicates), function(rep_idx) {
    run_one(f = f, seed_val = rep_idx * 17 + round(f * 1000))
  })
  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0) return(NULL)
  do.call(rbind, lapply(res_list, as.data.frame))
})
names(all_results) <- as.character(f_values)

# ---- Summarize ---------------------------------------------------------------

summary_rows <- lapply(names(all_results), function(f_str) {
  df <- all_results[[f_str]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  data.frame(
    resistance_fraction   = as.numeric(f_str),
    n_resistant_avg       = mean(df$n_resistant),
    mean_auroc_cell       = round(mean(df$auroc_cell, na.rm = TRUE), 3),
    sd_auroc_cell         = round(sd(df$auroc_cell,   na.rm = TRUE), 3),
    mean_auroc_clone      = round(mean(df$auroc_clone, na.rm = TRUE), 3),
    mean_correlation_Z    = round(mean(df$cor_z,        na.rm = TRUE), 3),
    mean_sensitivity      = round(mean(df$sensitivity,  na.rm = TRUE), 3),
    mean_specificity      = round(mean(df$specificity,  na.rm = TRUE), 3),
    mean_jaccard_features = round(mean(df$jaccard_features, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, Filter(Negate(is.null), summary_rows))

cat("\n=== Summary: CYFER performance across resistance prevalence ===\n")
print(summary_df)

# ---- Also run naive comparison (group cells by clone size) -------------------

run_naive <- function(f, seed_val) {
  set.seed(seed_val)
  n_resistant    <- max(1, round(f * n_cells))
  n_nonresistant <- n_cells - n_resistant
  is_resistant   <- c(rep(TRUE, n_resistant), rep(FALSE, n_nonresistant))
  cell_names     <- paste0("cell:", seq_len(n_cells))

  causal_mat <- rbind(
    matrix(rnorm(n_resistant    * n_causal, mu_resistant,  sigma_expr), nrow = n_resistant),
    matrix(rnorm(n_nonresistant * n_causal, mu_background, sigma_expr), nrow = n_nonresistant)
  )
  noise_mat  <- matrix(rnorm(n_cells * (n_features - n_causal), 0, 1), nrow = n_cells)
  X          <- cbind(causal_mat, noise_mat)
  rownames(X) <- cell_names
  colnames(X) <- feat_names
  X          <- scale(X)

  true_Z      <- ifelse(is_resistant, Z_resistant, Z_nonresistant) + rnorm(n_cells, 0, Z_noise)
  names(true_Z) <- cell_names
  cells_per_clone <- n_cells %/% n_clones
  clone_ids  <- rep(seq_len(n_clones), each = cells_per_clone)[seq_len(n_cells)]
  shuffle_idx <- sample(n_cells)
  clone_ids  <- clone_ids[order(shuffle_idx)]
  clone_labels <- paste0("clone:", clone_ids)

  cell_future  <- rpois(n_cells, lambda = exp(true_Z))
  names(cell_future) <- cell_names
  lineage_future_count <- tapply(cell_future, clone_labels, sum)

  # Naive score: assign cells the log-count of their clone
  naive_score  <- log(lineage_future_count[clone_labels] + 1)
  names(naive_score) <- cell_names

  auc_val <- auroc(naive_score, as.integer(is_resistant))
  list(f = f, auroc_naive = auc_val)
}

message("Running naive comparison...")
naive_results <- lapply(f_values, function(f) {
  res_list <- lapply(seq_len(n_replicates), function(rep_idx) {
    run_naive(f, seed_val = rep_idx * 17 + round(f * 1000))
  })
  mean(sapply(res_list, function(r) r$auroc_naive), na.rm = TRUE)
})

summary_df$mean_auroc_cell_naive <- round(unlist(naive_results), 3)

cat("\n=== Comparison: CYFER vs. naive (clone-level) cell AUROC ===\n")
print(summary_df[, c("resistance_fraction", "n_resistant_avg",
                     "mean_auroc_cell", "mean_auroc_cell_naive",
                     "mean_sensitivity", "mean_specificity")])

# ---- Save results ------------------------------------------------------------

saveRDS(
  list(all_results = all_results,
       summary     = summary_df,
       params      = list(n_cells = n_cells, n_clones = n_clones,
                          n_features = n_features, n_causal = n_causal,
                          f_values = f_values)),
  file = "sim2_rare_resistance_results.rds"
)

message("sim2 complete. Results saved to sim2_rare_resistance_results.rds")
