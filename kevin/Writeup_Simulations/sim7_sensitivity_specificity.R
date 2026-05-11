# ==============================================================================
# sim7_sensitivity_specificity.R
#
# GOAL: Provide a comprehensive sensitivity/specificity analysis comparing
# CYFER vs. the naive clone-size-based differential expression approach for
# identifying molecular features associated with fate potential, including
# formal statistical tests and ROC/precision-recall curves.
#
# REVIEWER CONCERNS:
#   Comment 7 — "The differential expression analysis in Figure 3E-G is not
#                clear: why is there no p-value associated with the CYFER
#                analysis? Additionally, what is the specificity and sensitivity
#                of these approaches in this DE analysis?"
#   Comment 1  — The simulation section should show "fits work when the model
#                is inaccurate" and provide quantitative metrics.
#
# DESIGN:
#   Based on the paper's existing priming and plasticity simulation framework,
#   but with added quantitative metrics and formal statistical comparisons.
#
#   Data generation (extending paper's simulation):
#     - n=1500 cells, p=50 features, K=60 clones.
#     - 10 features are truly causal (β ≠ 0); 40 are null (β = 0).
#     - Priming scenario: between-clone variance in fate potential is high.
#     - Plasticity scenario: within-clone variance in fate potential is high.
#
#   Methods compared:
#     (A) CYFER: estimate fate potential per cell; correlate features with
#         estimated fate potential; p-values from correlation test.
#     (B) Naive (DE by clone fate): divide clones into "high" and "low" groups
#         by observed t2 count; run Wilcoxon test for each feature.
#     (C) Clone-mean approach: compute mean fate potential per clone; test
#         correlation between clone mean expression and clone fate.
#
#   Metrics at each threshold:
#     - TP, FP, TN, FN (binary classification of each feature)
#     - Sensitivity (recall), Specificity, Precision, F1, Jaccard
#     - AUROC and AUPRC for ranking-based comparison
#     - p-value from permutation test (shuffling fate labels) to assess
#       whether the difference between CYFER and naive is statistically
#       significant.
#
#   This simulation is run for both priming and plasticity scenarios.
#
# EXPECTED RESULTS:
#   In the priming scenario, both methods should perform well because high-fate
#   cells are uniformly spread within each clone (so clone-level assignment
#   works). In the plasticity scenario, CYFER should substantially outperform
#   the naive method, because naive DE conflates the rare high-fate cells with
#   the majority of low-fate cells in the same clone.
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

jaccard <- function(a, b) {
  u <- union(a, b); if (length(u) == 0) return(0)
  length(intersect(a, b)) / length(u)
}

auroc <- function(scores, binary_labels) {
  n_pos <- sum(binary_labels == 1); n_neg <- sum(binary_labels == 0)
  if (n_pos == 0 || n_neg == 0) return(0.5)
  w <- wilcox.test(scores[binary_labels == 1], scores[binary_labels == 0],
                   alternative = "greater")
  as.numeric(w$statistic) / (n_pos * n_neg)
}

# Precision-recall AUC (trapezoid)
auprc <- function(scores, binary_labels) {
  thresholds <- sort(unique(scores), decreasing = TRUE)
  pts <- t(sapply(thresholds, function(t) {
    tp <- sum(scores >= t & binary_labels == 1)
    fp <- sum(scores >= t & binary_labels == 0)
    fn <- sum(scores <  t & binary_labels == 1)
    prec <- tp / max(tp + fp, 1)
    rec  <- tp / max(tp + fn, 1)
    c(rec = rec, prec = prec)
  }))
  # Trapezoid rule
  n <- nrow(pts)
  if (n < 2) return(pts[1, "prec"] * pts[1, "rec"])
  sum(diff(pts[, "rec"]) * (pts[-n, "prec"] + pts[-1, "prec"]) / 2)
}

compute_metrics_at_threshold <- function(p_values, true_causal, alpha = 0.05) {
  called <- names(p_values)[p_values < alpha]
  all_names <- names(p_values)
  tp <- length(intersect(called, true_causal))
  fp <- length(setdiff(called, true_causal))
  fn <- length(setdiff(true_causal, called))
  tn <- length(setdiff(all_names, union(called, true_causal)))
  sens <- tp / max(tp + fn, 1)
  spec <- tn / max(tn + fp, 1)
  prec <- tp / max(tp + fp, 1)
  f1   <- 2 * prec * sens / max(prec + sens, 1e-10)
  jac  <- jaccard(called, true_causal)
  c(TP = tp, FP = fp, TN = tn, FN = fn,
    Sensitivity = sens, Specificity = spec, Precision = prec, F1 = f1,
    Jaccard = jac, N_called = length(called))
}

# ---- Simulation parameters ---------------------------------------------------

n_cells    <- 1200
n_clones   <- 60
n_features <- 50
n_causal   <- 10

feat_names     <- paste0("f", seq_len(n_features))
causal_features <- feat_names[seq_len(n_causal)]
null_features   <- feat_names[(n_causal + 1):n_features]

# True coefficients
true_beta <- c(rep(1.5, n_causal), rep(0, n_features - n_causal))
names(true_beta) <- feat_names

sigma_between <- 1.0
sigma_within_priming   <- 0.3   # low within-clone variance → priming
sigma_within_plasticity <- 1.5  # high within-clone variance → plasticity

n_cells_per_clone <- n_cells %/% n_clones
clone_ids    <- rep(seq_len(n_clones), each = n_cells_per_clone)[seq_len(n_cells)]
clone_labels <- paste0("clone:", clone_ids)

n_replicates  <- 10
alpha_thresh  <- 0.05
n_permutations <- 100    # for permutation test

# ---- Methods to apply --------------------------------------------------------

# Method A: CYFER — correlate features with estimated fate potential
method_cyfer <- function(X, clone_sub, lfc_sub, true_Z) {
  n_folds_cv <- min(3, length(unique(clone_sub)) - 1)
  if (n_folds_cv < 2) return(NULL)
  tryCatch({
    fr <- cyfer(X, clone_sub, lfc_sub, lambda_initial = 1,
                lambda_sequence_length = 10, num_folds = n_folds_cv, verbose = 0)
    ff <- cyfer_finalize(X, clone_sub, fr, lfc_sub)
    Z_hat <- as.numeric(X %*% ff$coefficient_vec[-1]) + ff$coefficient_vec[1]
    # p-values: Pearson correlation test of each feature with Z_hat
    p_vals <- sapply(seq_len(ncol(X)), function(j) {
      if (sd(X[, j]) < 1e-10) return(1)
      ct <- cor.test(X[, j], Z_hat)
      ct$p.value
    })
    names(p_vals) <- colnames(X)
    cor_z <- cor(true_Z, Z_hat)
    list(p_values = p_vals, cor_z = cor_z, Z_hat = Z_hat)
  }, error = function(e) NULL)
}

# Method B: Naive (Wilcoxon by clone fate group: top 25% vs bottom 25% by count)
method_naive <- function(X, clone_sub, lfc_sub) {
  clone_score  <- lfc_sub[clone_sub]
  high_thresh  <- quantile(lfc_sub, 0.75)
  low_thresh   <- quantile(lfc_sub, 0.25)
  high_cells   <- which(clone_score >= high_thresh)
  low_cells    <- which(clone_score <= low_thresh)
  if (length(high_cells) < 2 || length(low_cells) < 2) return(NULL)
  p_vals <- sapply(seq_len(ncol(X)), function(j) {
    tryCatch(wilcox.test(X[high_cells, j], X[low_cells, j])$p.value,
             error = function(e) 1)
  })
  names(p_vals) <- colnames(X)
  list(p_values = p_vals)
}

# Method C: Clone-mean correlation (correlate each feature's clone mean with clone count)
method_clone_mean <- function(X, clone_sub, lfc_sub) {
  clone_means  <- t(sapply(unique(clone_sub), function(cl) {
    colMeans(X[clone_sub == cl, , drop = FALSE])
  }))
  rownames(clone_means) <- unique(clone_sub)
  lfc_vec <- lfc_sub[rownames(clone_means)]
  p_vals <- sapply(seq_len(ncol(clone_means)), function(j) {
    if (sd(clone_means[, j]) < 1e-10) return(1)
    cor.test(clone_means[, j], log(lfc_vec + 1))$p.value
  })
  names(p_vals) <- colnames(X)
  list(p_values = p_vals)
}

# ---- Run simulation for one scenario and replicate ---------------------------

run_one_replicate <- function(scenario, seed_val) {
  print(paste0("Seed: ", seed_val))
  set.seed(seed_val)

  sigma_within <- if (scenario == "priming") sigma_within_priming else sigma_within_plasticity

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
  lfc <- tapply(cell_future, clone_labels, sum)

  valid_clones <- names(lfc[lfc > 0]); if (length(valid_clones) < 5) return(NULL)
  keep    <- which(clone_labels %in% valid_clones)
  X_sub   <- X[keep, , drop = FALSE]
  Z_sub   <- true_Z[keep]
  cl_sub  <- clone_labels[keep]
  lfc_sub <- lfc[valid_clones]

  cs <- table(cl_sub); valid2 <- names(cs[cs >= 2])
  if (length(valid2) < 5) return(NULL)
  k2 <- which(cl_sub %in% valid2)
  X_sub <- X_sub[k2, , drop=FALSE]; Z_sub <- Z_sub[k2]; cl_sub <- cl_sub[k2]
  lfc_sub <- lfc_sub[valid2]

  # Run all three methods
  res_cyfer      <- method_cyfer(X_sub, cl_sub, lfc_sub, Z_sub)
  res_naive      <- method_naive(X_sub, cl_sub, lfc_sub)
  res_clone_mean <- method_clone_mean(X_sub, cl_sub, lfc_sub)

  if (is.null(res_cyfer) || is.null(res_naive)) return(NULL)

  # Metrics at alpha = 0.05
  m_cyfer      <- compute_metrics_at_threshold(res_cyfer$p_values, causal_features)
  m_naive      <- compute_metrics_at_threshold(res_naive$p_values, causal_features)
  m_clone_mean <- if (!is.null(res_clone_mean))
    compute_metrics_at_threshold(res_clone_mean$p_values, causal_features) else NA

  # AUROC (ranking-based)
  true_binary  <- as.integer(feat_names %in% causal_features)
  scores_cyfer <- -log10(res_cyfer$p_values + 1e-300)
  scores_naive <- -log10(res_naive$p_values + 1e-300)
  auc_cyfer    <- auroc(scores_cyfer, true_binary)
  auc_naive    <- auroc(scores_naive, true_binary)
  prc_cyfer    <- auprc(scores_cyfer, true_binary)
  prc_naive    <- auprc(scores_naive, true_binary)

  list(
    scenario     = scenario,
    seed         = seed_val,
    cor_z_cyfer  = if (!is.null(res_cyfer)) res_cyfer$cor_z else NA,
    metrics_cyfer      = m_cyfer,
    metrics_naive      = m_naive,
    metrics_clone_mean = m_clone_mean,
    auroc_cyfer  = auc_cyfer,
    auroc_naive  = auc_naive,
    auprc_cyfer  = prc_cyfer,
    auprc_naive  = prc_naive
  )
}

# ---- Run both scenarios ------------------------------------------------------

scenarios <- c("priming", "plasticity")

message("Running ", n_replicates, " replicates x ", length(scenarios), " scenarios...")

all_results <- lapply(scenarios, function(sc) {
  message("  Scenario: ", sc)
  reps <- lapply(seq_len(n_replicates), function(rep_idx) {
    run_one_replicate(sc, seed_val = rep_idx * 37 + ifelse(sc == "priming", 0, 1000))
  })
  Filter(Negate(is.null), reps)
})
names(all_results) <- scenarios

# ---- Summarize metrics -------------------------------------------------------

summarize_scenario <- function(sc_name, reps) {
  extract <- function(reps, field) sapply(reps, function(r) r[[field]])
  extract_metric <- function(reps, method, metric) {
    sapply(reps, function(r) r[[method]][[metric]])
  }

  data.frame(
    scenario          = sc_name,
    n_converged       = length(reps),
    cor_z_cyfer       = round(mean(extract(reps, "cor_z_cyfer"),  na.rm = TRUE), 3),
    auroc_cyfer       = round(mean(extract(reps, "auroc_cyfer"),  na.rm = TRUE), 3),
    auroc_naive       = round(mean(extract(reps, "auroc_naive"),  na.rm = TRUE), 3),
    auprc_cyfer       = round(mean(extract(reps, "auprc_cyfer"),  na.rm = TRUE), 3),
    auprc_naive       = round(mean(extract(reps, "auprc_naive"),  na.rm = TRUE), 3),
    sens_cyfer        = round(mean(extract_metric(reps, "metrics_cyfer", "Sensitivity"), na.rm=TRUE), 3),
    spec_cyfer        = round(mean(extract_metric(reps, "metrics_cyfer", "Specificity"), na.rm=TRUE), 3),
    prec_cyfer        = round(mean(extract_metric(reps, "metrics_cyfer", "Precision"),   na.rm=TRUE), 3),
    f1_cyfer          = round(mean(extract_metric(reps, "metrics_cyfer", "F1"),          na.rm=TRUE), 3),
    jac_cyfer         = round(mean(extract_metric(reps, "metrics_cyfer", "Jaccard"),     na.rm=TRUE), 3),
    sens_naive        = round(mean(extract_metric(reps, "metrics_naive", "Sensitivity"), na.rm=TRUE), 3),
    spec_naive        = round(mean(extract_metric(reps, "metrics_naive", "Specificity"), na.rm=TRUE), 3),
    prec_naive        = round(mean(extract_metric(reps, "metrics_naive", "Precision"),   na.rm=TRUE), 3),
    f1_naive          = round(mean(extract_metric(reps, "metrics_naive", "F1"),          na.rm=TRUE), 3),
    jac_naive         = round(mean(extract_metric(reps, "metrics_naive", "Jaccard"),     na.rm=TRUE), 3),
    stringsAsFactors  = FALSE
  )
}

summary_list <- lapply(scenarios, function(sc) {
  summarize_scenario(sc, all_results[[sc]])
})
summary_df <- do.call(rbind, summary_list)

cat("\n=== Main results table ===\n")
cat("AUROC (area under ROC curve), AUPRC (area under precision-recall curve)\n\n")
print(summary_df[, c("scenario", "n_converged", "cor_z_cyfer",
                     "auroc_cyfer", "auroc_naive",
                     "sens_cyfer",  "spec_cyfer",  "jac_cyfer",
                     "sens_naive",  "spec_naive",  "jac_naive")])

# ---- Permutation test: Is CYFER significantly better than naive? -------------

message("\nRunning permutation tests...")

permutation_test <- lapply(scenarios, function(sc) {
  reps <- all_results[[sc]]

  # Observed difference in AUROC
  obs_auc_cyfer <- sapply(reps, `[[`, "auroc_cyfer")
  obs_auc_naive <- sapply(reps, `[[`, "auroc_naive")
  obs_diff      <- mean(obs_auc_cyfer - obs_auc_naive, na.rm = TRUE)

  # Permutation: shuffle which AUC came from which method
  perm_diffs <- replicate(n_permutations, {
    combined <- c(obs_auc_cyfer, obs_auc_naive)
    n        <- length(obs_auc_cyfer)
    shuffled <- sample(combined, length(combined))
    mean(shuffled[seq_len(n)] - shuffled[(n+1):(2*n)], na.rm = TRUE)
  })
  p_value <- mean(perm_diffs >= obs_diff)

  data.frame(scenario = sc, obs_diff_auroc = round(obs_diff, 4),
             p_value = round(p_value, 4), stringsAsFactors = FALSE)
})
perm_df <- do.call(rbind, permutation_test)

cat("\n=== Permutation test: CYFER AUROC vs naive AUROC ===\n")
print(perm_df)

# ---- P-value calibration check (false positive rate under null) --------------

message("\nRunning null calibration (checking false positive rates)...")

# Under the null: all features have beta = 0
null_run <- function(seed_val) {
  print(paste0("Seed: ", seed_val))
  set.seed(seed_val)
  X_null <- matrix(rnorm(n_cells * n_features), nrow = n_cells)
  rownames(X_null) <- paste0("cell:", seq_len(n_cells))
  colnames(X_null) <- feat_names
  X_null <- scale(X_null)

  true_Z_null <- rnorm(n_cells); names(true_Z_null) <- rownames(X_null)
  cell_future <- rpois(n_cells, exp(true_Z_null)); names(cell_future) <- rownames(X_null)
  lfc <- tapply(cell_future, clone_labels, sum)

  valid <- names(lfc[lfc > 0]); if (length(valid) < 5) return(NULL)
  keep  <- which(clone_labels %in% valid)
  X_s   <- X_null[keep, , drop = FALSE]; cl_s <- clone_labels[keep]; lfc_s <- lfc[valid]

  cs <- table(cl_s); valid2 <- names(cs[cs >= 2])
  if (length(valid2) < 5) return(NULL)
  k2 <- which(cl_s %in% valid2)
  X_s <- X_s[k2, , drop=FALSE]; cl_s <- cl_s[k2]; lfc_s <- lfc_s[valid2]

  nf <- min(3, length(unique(cl_s)) - 1); if (nf < 2) return(NULL)
  tryCatch({
    fr <- cyfer(X_s, cl_s, lfc_s, lambda_initial = 1,
                lambda_sequence_length = 10, num_folds = nf, verbose = 0)
    ff <- cyfer_finalize(X_s, cl_s, fr, lfc_s)
    Z_hat <- as.numeric(X_s %*% ff$coefficient_vec[-1]) + ff$coefficient_vec[1]
    p_vals <- sapply(seq_len(ncol(X_s)), function(j) cor.test(X_s[, j], Z_hat)$p.value)
    mean(p_vals < alpha_thresh)   # false positive rate
  }, error = function(e) NA)
}

null_fpr <- sapply(seq_len(10), null_run)
cat(sprintf("\nNull calibration: expected FPR=%.2f, observed mean FPR=%.3f (SD=%.3f)\n",
            alpha_thresh, mean(null_fpr, na.rm = TRUE), sd(null_fpr, na.rm = TRUE)))

# ---- Save all results --------------------------------------------------------

filepath <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup_Simulations/"
saveRDS(
  list(summary        = summary_df,
       permutation    = perm_df,
       all_results    = all_results,
       null_fpr       = null_fpr,
       params = list(n_cells = n_cells, n_clones = n_clones,
                     n_features = n_features, n_causal = n_causal,
                     alpha = alpha_thresh, n_replicates = n_replicates)),
  file = paste0(filepath, "sim7_sensitivity_specificity_results.rds")
)

message("sim7 complete. Results saved to sim7_sensitivity_specificity_results.rds")
