# ==============================================================================
# sim3_heritability.R
#
# GOAL: Test CYFER's performance across varying levels of phenotypic heritability —
# i.e., the degree to which a cell's molecular features (expression) predict its
# fate potential, and the degree to which sister cells within a clone share
# similar expression profiles.
#
# REVIEWER CONCERNS:
#   Comment 1c — "realistically modelling phenotype heritability."
#   Comment 2  — The 'plasticity' scenario could equally be explained by mutation
#                within a clone after barcoding (not traditional plasticity).
#   Comment 3  — "generating synthetic data to preserve phenotype heritability
#                across lineages would be better validation of CYFER."
#
# DESIGN:
#   Two axes of heritability are simulated independently and jointly:
#
#   (1) Feature heritability (h2_feature): How strongly do shared clone-level
#       gene expression programs predict fate?
#       - Controlled by the ratio of between-clone to within-clone expression
#         variance. h2_feature = sigma_between^2 / (sigma_between^2 + sigma_within^2)
#       - High h2_feature: sister cells in the same clone have similar expression
#         → gene expression is a reliable lineage marker.
#       - Low h2_feature: sister cells diverge in expression → expression is
#         largely noise relative to lineage.
#
#   (2) Fate heritability (h2_fate): How much of the variance in fate potential
#       is explained by the shared clone mean (vs. cell-level stochasticity)?
#       - Controlled by residual fate noise (sigma_fate).
#       - High h2_fate: cell fate is tightly linked to its expression.
#       - Low h2_fate: fate is mostly random (epistatic / environmental).
#
#   Simulations span a grid of (h2_feature, h2_fate) values.
#   CYFER is applied to each, and performance (correlation, Jaccard) is reported.
#
# EXPECTED RESULTS:
#   CYFER performance should increase with higher h2_feature (clearer expression
#   signal for clones) and higher h2_fate (fate more predictable from expression).
#   Even at moderate heritability, CYFER should outperform naive clone-level DE,
#   because it explicitly weights each cell by its fate potential.
# ==============================================================================

rm(list=ls())
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

top_features_by_cor <- function(feature_mat, target, k) {
  cors <- apply(feature_mat, 2, function(col) {
    if (sd(col) < 1e-10) return(0)
    cor(col, target)
  })
  names(sort(abs(cors), decreasing = TRUE))[seq_len(min(k, ncol(feature_mat)))]
}

# ---- Simulation parameters ---------------------------------------------------

n_cells      <- 1200
n_clones     <- 50
n_features   <- 25
n_causal     <- 5       # features that have non-zero coefficient in true model

# True coefficients (sparse)
true_beta      <- c(rep(1.5, n_causal), rep(0, n_features - n_causal))
feat_names     <- paste0("f", seq_len(n_features))
names(true_beta) <- feat_names
causal_features <- feat_names[seq_len(n_causal)]

# Cell-to-clone assignment
cells_per_clone <- n_cells %/% n_clones
clone_ids       <- rep(seq_len(n_clones), each = cells_per_clone)[seq_len(n_cells)]
clone_labels    <- paste0("clone:", clone_ids)

# Heritability grid
# h2_feature controls how heritable expression is within a clone
# h2_fate controls how much fate is explained by expression
h2_feature_values <- c(0.05, 0.20, 0.50, 0.80, 0.95)
h2_fate_values    <- c(0.10, 0.30, 0.60, 0.90)

n_replicates <- 10

# ---- Core simulation function ------------------------------------------------

# Convert h2 to sigma ratio: h2 = sigma_b^2 / (sigma_b^2 + sigma_w^2)
# If we fix sigma_b = 1, then sigma_w = sqrt((1 - h2) / h2)
h2_to_sigma_within <- function(h2, sigma_between = 1.0) {
  sigma_between * sqrt((1 - h2) / max(h2, 1e-6))
}

run_heritability_sim <- function(h2_feature, h2_fate, seed_val) {
  print(paste0("Seed: ", seed_val))
  set.seed(seed_val)

  sigma_between <- 1.0
  sigma_within  <- h2_to_sigma_within(h2_feature, sigma_between)

  # Generate hierarchical expression matrix
  clone_centers <- mvrnorm(n_clones,
                           mu    = rep(0, n_features),
                           Sigma = sigma_between^2 * diag(n_features))
  rownames(clone_centers) <- paste0("clone:", seq_len(n_clones))
  colnames(clone_centers) <- feat_names

  X <- t(sapply(seq_len(n_cells), function(i) {
    k <- clone_ids[i]
    clone_centers[k, ] +
      mvrnorm(1, mu = rep(0, n_features),
              Sigma = sigma_within^2 * diag(n_features))
  }))
  rownames(X) <- paste0("cell:", seq_len(n_cells))
  colnames(X) <- feat_names
  X <- scale(X)

  # True fate potential: linear in expression + cell-level noise
  # h2_fate = Var(beta^T X_mean) / (Var(beta^T X_mean) + sigma_fate^2)
  # We compute sigma_fate from h2_fate
  Z_from_expr <- as.numeric(X %*% true_beta)
  var_z_expr  <- var(Z_from_expr)
  sigma_fate  <- sqrt(var_z_expr * (1 - h2_fate) / max(h2_fate, 1e-6))
  true_Z      <- Z_from_expr + rnorm(n_cells, 0, sigma_fate)
  names(true_Z) <- rownames(X)
  true_Z         <- true_Z - mean(true_Z)   # center

  # Generate observed future clone sizes
  cell_future <- rpois(n_cells, lambda = exp(true_Z))
  names(cell_future) <- rownames(X)
  lineage_future_count <- tapply(cell_future, clone_labels, sum)

  valid_clones <- names(lineage_future_count[lineage_future_count > 0])
  if (length(valid_clones) < 5) return(NULL)

  keep_idx  <- which(clone_labels %in% valid_clones)
  X_sub     <- X[keep_idx, , drop = FALSE]
  Z_sub     <- true_Z[keep_idx]
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

    cor_z     <- cor(Z_sub, Z_hat)
    est_feat  <- top_features_by_cor(X_sub, Z_hat, k = n_causal)
    jac_val   <- jaccard(causal_features, est_feat)
    gini_val  <- gini_coef(exp(Z_hat))

    list(h2_feature = h2_feature, h2_fate = h2_fate, seed = seed_val,
         cor_z = cor_z, jaccard = jac_val, gini = gini_val,
         n_valid_clones = length(valid_clones))
  }, error = function(e) NULL)
}

# ---- Naive baseline (group cells by clone future count) ----------------------

run_naive_heritability <- function(h2_feature, h2_fate, seed_val) {
  print(paste0("Seed: ", seed_val))
  set.seed(seed_val)

  sigma_between <- 1.0
  sigma_within  <- h2_to_sigma_within(h2_feature, sigma_between)

  clone_centers <- mvrnorm(n_clones, mu = rep(0, n_features),
                           Sigma = sigma_between^2 * diag(n_features))
  colnames(clone_centers) <- feat_names

  X <- t(sapply(seq_len(n_cells), function(i) {
    k <- clone_ids[i]
    clone_centers[k, ] +
      mvrnorm(1, mu = rep(0, n_features), Sigma = sigma_within^2 * diag(n_features))
  }))
  colnames(X) <- feat_names
  X <- scale(X)

  Z_from_expr <- as.numeric(X %*% true_beta)
  var_z_expr  <- var(Z_from_expr)
  sigma_fate  <- sqrt(var_z_expr * (1 - h2_fate) / max(h2_fate, 1e-6))
  true_Z      <- Z_from_expr + rnorm(n_cells, 0, sigma_fate) - mean(Z_from_expr)

  cell_future <- rpois(n_cells, lambda = exp(true_Z))
  lineage_future_count <- tapply(cell_future, clone_labels, sum)

  # Naive score: cells in large clones get high score
  clone_score  <- log(lineage_future_count[clone_labels] + 1)
  cor_z_naive  <- cor(true_Z, clone_score)

  est_feat_naive <- top_features_by_cor(X, clone_score, k = n_causal)
  jac_naive      <- jaccard(causal_features, est_feat_naive)

  list(h2_feature = h2_feature, h2_fate = h2_fate,
       cor_z_naive = cor_z_naive, jaccard_naive = jac_naive)
}

# ---- Run all combinations ----------------------------------------------------

message("Running heritability grid (", length(h2_feature_values), " x ",
        length(h2_fate_values), " x ", n_replicates, " replicates)...")

grid <- expand.grid(h2_feature = h2_feature_values,
                    h2_fate    = h2_fate_values,
                    stringsAsFactors = FALSE)

all_cyfer <- lapply(seq_len(nrow(grid)), function(gi) {
  h2f   <- grid$h2_feature[gi]
  h2fa  <- grid$h2_fate[gi]
  message("  h2_feature=", h2f, ", h2_fate=", h2fa)

  cyfer_res <- lapply(seq_len(n_replicates), function(rep_idx) {
    run_heritability_sim(h2f, h2fa, seed_val = rep_idx * 31 + gi * 7)
  })
  cyfer_res <- Filter(Negate(is.null), cyfer_res)

  naive_res <- lapply(seq_len(n_replicates), function(rep_idx) {
    run_naive_heritability(h2f, h2fa, seed_val = rep_idx * 31 + gi * 7)
  })
  naive_res <- Filter(Negate(is.null), naive_res)

  data.frame(
    h2_feature        = h2f,
    h2_fate           = h2fa,
    mean_cor_cyfer    = round(mean(sapply(cyfer_res, `[[`, "cor_z"),    na.rm = TRUE), 3),
    sd_cor_cyfer      = round(sd(  sapply(cyfer_res, `[[`, "cor_z"),    na.rm = TRUE), 3),
    mean_jac_cyfer    = round(mean(sapply(cyfer_res, `[[`, "jaccard"),  na.rm = TRUE), 3),
    mean_cor_naive    = round(mean(sapply(naive_res,  `[[`, "cor_z_naive"),   na.rm = TRUE), 3),
    mean_jac_naive    = round(mean(sapply(naive_res,  `[[`, "jaccard_naive"), na.rm = TRUE), 3),
    stringsAsFactors  = FALSE
  )
})

summary_df <- do.call(rbind, all_cyfer)
summary_df$cor_improvement <- summary_df$mean_cor_cyfer - summary_df$mean_cor_naive
summary_df$jac_improvement <- summary_df$mean_jac_cyfer - summary_df$mean_jac_naive

cat("\n=== Heritability grid results ===\n")
print(summary_df)

cat("\n=== CYFER advantage (cor improvement over naive) ===\n")
# Show as matrix: rows = h2_feature, cols = h2_fate
for (h2f in h2_feature_values) {
  row_df <- summary_df[summary_df$h2_feature == h2f, ]
  cat(sprintf("h2_feature=%.2f | CYFER cor: %s | Naive cor: %s\n",
              h2f,
              paste(sprintf("%.3f", row_df$mean_cor_cyfer), collapse=" "),
              paste(sprintf("%.3f", row_df$mean_cor_naive),  collapse=" ")))
}

# ---- Save results ------------------------------------------------------------

filepath <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup_Simulations/"
saveRDS(
  list(summary = summary_df,
       grid    = grid,
       params  = list(n_cells = n_cells, n_clones = n_clones,
                      n_features = n_features, n_causal = n_causal,
                      h2_feature_values = h2_feature_values,
                      h2_fate_values    = h2_fate_values)),
  file = paste0(filepath, "sim3_heritability_results.rds")
)

message("sim3 complete. Results saved to sim3_heritability_results.rds")
