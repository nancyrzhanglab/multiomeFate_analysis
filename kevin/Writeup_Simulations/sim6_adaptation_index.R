# ==============================================================================
# sim6_adaptation_index.R
#
# GOAL: Test whether CYFER's fate-potential weighting in the adaptation index
# is robust to confounding from "dying" cells that show stress-like expression
# signatures but have negligible long-term fate potential.
#
# REVIEWER CONCERN: Comment 8 — "one would still expect the inverse correlation
# between selection and adaptation; consider this counter-example: Suppose that
# a population contains 100 cells, with 1% 'primed' for resistance and the
# other 99% dying off... at 5 weeks, the population would present a low
# 'adaptation' signature as just the original 1% of primed cells would remain.
# Indeed, the authors present data on p.10 line 389 that the adapting cells
# are associated with stress-like gene signatures. Can the authors offer any
# evidence that this temporal bias and actively dying cells are not driving
# the adaptation signal in the dataset?"
#
# DESIGN:
#   We simulate a clone that contains two cell types:
#     (1) "Surviving" cells (fraction f_survive): high fate potential Z_high,
#         stable or slightly changed expression from t1 to t2.
#     (2) "Dying" cells (fraction 1 - f_survive): low/negative fate potential
#         Z_low, but large expression change from t1 to t2 (stress-like).
#         These cells would dominate an UNWEIGHTED adaptation index.
#
#   The adaptation index is defined as:
#     I_adapt = (1 / sum_i(w_i)) * sum_i(w_i * d_i)
#   where w_i is a weighting (fate potential in CYFER, or uniform in naive),
#   and d_i is the distance in expression space from cell i to a "reference"
#   t2 population center.
#
#   We compare:
#     (A) CYFER-weighted adaptation index (w_i = estimated fate potential)
#     (B) Naive unweighted adaptation index (w_i = 1 for all cells)
#     (C) Oracle weighted adaptation index (w_i = true fate potential)
#
#   Key variable: the fraction of dying cells, varied from 10% to 99%.
#   Key metric: does the adaptation index reflect the surviving cells' state
#               change, not the dying cells' stress signature?
#
# EXPECTED RESULTS:
#   The naive unweighted adaptation index will be inflated by dying cells
#   (which have large d_i due to stress-induced expression changes). The
#   CYFER-weighted index should correctly downweight dying cells, producing
#   a lower, more accurate adaptation score that reflects only the state
#   changes of cells that actually contribute to future growth.
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

# Weighted adaptation index: weighted mean of distances d_i
adaptation_index <- function(d_vec, weights) {
  weights <- pmax(weights, 0)
  if (sum(weights) < 1e-10) return(mean(d_vec))
  sum(weights * d_vec) / sum(weights)
}

fit_cyfer_simple <- function(X, clone_labels, lfc) {
  valid <- names(lfc[lfc > 0])
  if (length(valid) < 3) return(NULL)
  keep <- which(clone_labels %in% valid)
  X_s  <- X[keep, , drop = FALSE]; cl_s <- clone_labels[keep]; lfc_s <- lfc[valid]
  cs <- table(cl_s); valid2 <- names(cs[cs >= 2])
  if (length(valid2) < 3) return(NULL)
  k2 <- which(cl_s %in% valid2); X_s <- X_s[k2, , drop=FALSE]
  cl_s <- cl_s[k2]; lfc_s <- lfc_s[valid2]
  nf <- min(3, length(unique(cl_s)) - 1); if (nf < 2) return(NULL)
  tryCatch({
    fr <- cyfer(X_s, cl_s, lfc_s, lambda_initial = 1,
                lambda_sequence_length = 10, num_folds = nf, verbose = 0)
    cyfer_finalize(X_s, cl_s, fr, lfc_s)
  }, error = function(e) NULL)
}

# ---- Simulation parameters ---------------------------------------------------

n_cells_per_clone <- 200   # cells per clone at t1
n_clones          <- 60    # number of clones

n_features  <- 20          # expression dimensions
n_features_stress <- 5     # features that change in dying cells (stress signature)

# Fate potential parameters
Z_high <- 2.5    # log expected progeny for surviving cells
Z_low  <- -3.0   # log expected progeny for dying cells

# Expression parameters at t1
mu_survive_t1 <- c(rep(1.5, 5), rep(0, n_features - 5))   # survival signature
mu_dying_t1   <- c(rep(0, 5),   rep(0, n_features - 5))   # no distinctive t1 signature
sigma_expr    <- 0.5

# Expression change from t1 to t2 for dying cells (stress-induced)
delta_stress  <- c(rep(0, 5), rep(3.0, n_features_stress), rep(0, n_features - 5 - n_features_stress))
# Surviving cells: minimal expression change
delta_survive <- c(rep(0.5, 5), rep(0, n_features - 5))

# Fractions of dying cells to test
f_dying_values <- c(0.10, 0.25, 0.50, 0.75, 0.90, 0.99)

n_replicates <- 20

feat_names <- paste0("f", seq_len(n_features))

# ---- Core simulation for a single clone proportions -------------------------

run_one <- function(f_dying, seed_val) {
  set.seed(seed_val)

  n_dying   <- round(f_dying * n_cells_per_clone)
  n_survive <- n_cells_per_clone - n_dying
  is_dying  <- c(rep(TRUE, n_dying), rep(FALSE, n_survive))

  # Cell IDs and random clone assignment
  n_total      <- n_cells_per_clone * n_clones
  n_survive_tot <- n_survive * n_clones
  n_dying_tot  <- n_dying * n_clones

  clone_ids    <- rep(seq_len(n_clones), each = n_cells_per_clone)
  cell_types   <- rep(is_dying, times = n_clones)   # which cells are dying
  cell_names   <- paste0("cell:", seq_len(n_total))
  clone_labels <- paste0("clone:", clone_ids)

  # Expression at t1 (early time point)
  X_t1 <- t(sapply(seq_len(n_total), function(i) {
    if (cell_types[i]) {
      mu_dying_t1 + rnorm(n_features, 0, sigma_expr)
    } else {
      mu_survive_t1 + rnorm(n_features, 0, sigma_expr)
    }
  }))
  rownames(X_t1) <- cell_names; colnames(X_t1) <- feat_names
  X_t1 <- scale(X_t1)

  # Expression at t2 (late time point, simulated as a "target state")
  # Each cell's expected t2 state is its t1 state + type-specific shift
  X_t2_expected <- t(sapply(seq_len(n_total), function(i) {
    if (cell_types[i]) {
      X_t1[i, ] + delta_stress + rnorm(n_features, 0, 0.3)
    } else {
      X_t1[i, ] + delta_survive + rnorm(n_features, 0, 0.3)
    }
  }))
  rownames(X_t2_expected) <- cell_names; colnames(X_t2_expected) <- feat_names

  # Per-cell "distance" to t2 state (what adaptation index measures)
  d_vec <- sqrt(rowSums((X_t1 - X_t2_expected)^2))
  names(d_vec) <- cell_names

  # True fate potential
  true_Z <- ifelse(cell_types, Z_low, Z_high) + rnorm(n_total, 0, 0.3)
  names(true_Z) <- cell_names

  # Observed clone sizes at t2 (from Poisson model)
  cell_future <- rpois(n_total, pmax(exp(true_Z), 0.01))
  names(cell_future) <- cell_names
  lfc <- tapply(cell_future, clone_labels, sum)

  # Fit CYFER on t1 expression to estimate fate potential
  final_fit <- fit_cyfer_simple(X_t1, clone_labels, lfc)
  if (is.null(final_fit)) return(NULL)

  # Estimated fate potential (linear predictor from CYFER)
  Z_hat <- as.numeric(X_t1 %*% final_fit$coefficient_vec[-1]) +
           final_fit$coefficient_vec[1]
  names(Z_hat) <- cell_names
  exp_Z_hat <- exp(Z_hat)

  # ---- Adaptation indices ----

  # (A) CYFER-weighted: use estimated fate potential as weights
  # Compute per-clone, then average over clones weighted by clone size
  adapt_cyfer  <- tapply(seq_len(n_total), clone_labels, function(idx) {
    adaptation_index(d_vec[idx], exp_Z_hat[idx])
  })

  # (B) Naive unweighted: uniform weights
  adapt_naive  <- tapply(seq_len(n_total), clone_labels, function(idx) {
    mean(d_vec[idx])
  })

  # (C) Oracle: use true fate potential as weights
  adapt_oracle <- tapply(seq_len(n_total), clone_labels, function(idx) {
    adaptation_index(d_vec[idx], exp(true_Z[idx]))
  })

  # "True" adaptation should reflect surviving cells' state change, not dying cells'
  # Ground truth: adaptation = mean d_i for surviving cells only (oracle reference)
  adapt_survive_only <- tapply(seq_len(n_total), clone_labels, function(idx) {
    survive_idx <- idx[!cell_types[idx]]
    if (length(survive_idx) == 0) return(NA)
    mean(d_vec[survive_idx])
  })

  # Summary across clones (mean ± SD)
  bias_cyfer  <- mean(adapt_cyfer  - adapt_survive_only, na.rm = TRUE)
  bias_naive  <- mean(adapt_naive  - adapt_survive_only, na.rm = TRUE)
  bias_oracle <- mean(adapt_oracle - adapt_survive_only, na.rm = TRUE)

  # Correlation of adaptation index with true survive-only adaptation
  cor_cyfer  <- cor(adapt_cyfer,        adapt_survive_only, use = "complete.obs")
  cor_naive  <- cor(adapt_naive,         adapt_survive_only, use = "complete.obs")
  cor_oracle <- cor(adapt_oracle,        adapt_survive_only, use = "complete.obs")

  # Correlation with fate potential (confirms CYFER correctly identified Z)
  cor_fate <- cor(true_Z, Z_hat)

  list(
    f_dying      = f_dying,
    bias_cyfer   = bias_cyfer,  bias_naive = bias_naive,  bias_oracle = bias_oracle,
    cor_cyfer    = cor_cyfer,   cor_naive  = cor_naive,   cor_oracle  = cor_oracle,
    cor_fate     = cor_fate,
    mean_adapt_cyfer  = mean(adapt_cyfer,        na.rm = TRUE),
    mean_adapt_naive  = mean(adapt_naive,         na.rm = TRUE),
    mean_adapt_oracle = mean(adapt_oracle,        na.rm = TRUE),
    mean_adapt_true   = mean(adapt_survive_only,  na.rm = TRUE)
  )
}

# ---- Run all settings --------------------------------------------------------

message("Running adaptation index simulation across f_dying values...")

all_results <- lapply(f_dying_values, function(f_d) {
  message("  f_dying = ", f_d)
  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    run_one(f_d, seed_val = rep_idx * 29 + round(f_d * 100))
  })
  Filter(Negate(is.null), reps)
})
names(all_results) <- as.character(f_dying_values)

# ---- Summarize results -------------------------------------------------------

summary_rows <- lapply(names(all_results), function(f_str) {
  reps <- all_results[[f_str]]
  if (length(reps) == 0) return(NULL)
  data.frame(
    f_dying              = as.numeric(f_str),
    f_survive            = 1 - as.numeric(f_str),
    mean_bias_naive      = round(mean(sapply(reps, `[[`, "bias_naive"),  na.rm = TRUE), 4),
    mean_bias_cyfer      = round(mean(sapply(reps, `[[`, "bias_cyfer"),  na.rm = TRUE), 4),
    mean_bias_oracle     = round(mean(sapply(reps, `[[`, "bias_oracle"), na.rm = TRUE), 4),
    mean_cor_naive       = round(mean(sapply(reps, `[[`, "cor_naive"),   na.rm = TRUE), 3),
    mean_cor_cyfer       = round(mean(sapply(reps, `[[`, "cor_cyfer"),   na.rm = TRUE), 3),
    mean_cor_oracle      = round(mean(sapply(reps, `[[`, "cor_oracle"),  na.rm = TRUE), 3),
    mean_cor_fate        = round(mean(sapply(reps, `[[`, "cor_fate"),    na.rm = TRUE), 3),
    mean_adapt_naive_vs_truth = round(
      mean(sapply(reps, `[[`, "mean_adapt_naive") - sapply(reps, `[[`, "mean_adapt_true"),
           na.rm = TRUE), 4),
    mean_adapt_cyfer_vs_truth = round(
      mean(sapply(reps, `[[`, "mean_adapt_cyfer") - sapply(reps, `[[`, "mean_adapt_true"),
           na.rm = TRUE), 4),
    n_converged          = length(reps),
    stringsAsFactors     = FALSE
  )
})
summary_df <- do.call(rbind, Filter(Negate(is.null), summary_rows))

cat("\n=== Adaptation index: CYFER vs naive vs oracle ===\n")
cat("  bias = mean(estimated adaptation - true surviving-cell adaptation)\n")
cat("  Positive bias = overestimation (dying cells inflate the index)\n\n")
print(summary_df[, c("f_dying", "mean_bias_naive", "mean_bias_cyfer",
                     "mean_bias_oracle", "mean_cor_naive", "mean_cor_cyfer")])

cat("\n=== Adaptation index values by method ===\n")
print(summary_df[, c("f_dying", "mean_adapt_naive_vs_truth",
                     "mean_adapt_cyfer_vs_truth")])

# ---- Additional: Inflection analysis -----------------------------------------
# At what dying cell fraction does the naive approach break down badly?
cat("\n=== Naive vs CYFER adaptation bias ratio ===\n")
ratio_df <- data.frame(
  f_dying        = summary_df$f_dying,
  bias_ratio     = round(summary_df$mean_bias_naive / (summary_df$mean_bias_cyfer + 1e-6), 2)
)
print(ratio_df)

# ---- Save results ------------------------------------------------------------

saveRDS(
  list(all_results = all_results,
       summary     = summary_df,
       params = list(n_cells_per_clone = n_cells_per_clone,
                     n_clones = n_clones, n_features = n_features,
                     Z_high = Z_high, Z_low = Z_low,
                     f_dying_values = f_dying_values)),
  file = "sim6_adaptation_index_results.rds"
)

message("sim6 complete. Results saved to sim6_adaptation_index_results.rds")
