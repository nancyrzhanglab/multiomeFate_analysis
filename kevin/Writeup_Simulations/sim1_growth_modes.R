# ==============================================================================
# sim1_growth_modes.R
#
# GOAL: Assess CYFER's robustness when the true underlying growth model deviates
# from the assumed exponential (Poisson-log) model.
#
# REVIEWER CONCERN: Comment 1a — "model (a) different growth modes, linear,
# exponential, others, and also fit the 'wrong' growth model to each scenario
# to show the fits work when the model is inaccurate."
#
# DESIGN:
#   Synthetic embedding: n=1500 cells, p=30 features, K=60 clones.
#   True fate potential: Z_i = beta^T X_i (linear combination of features).
#   Four growth models generate observed clone sizes from Z_i:
#     (A) Exponential: count_i ~ Poisson(exp(Z_i))          [CYFER's model]
#     (B) Linear:      count_i ~ Poisson(max(a + b*Z_i, 0)) [linear scaling]
#     (C) Logistic:    count_i ~ Poisson(C / (1+exp(-Z_i))) [saturating/bounded]
#     (D) Power-law:   count_i ~ Poisson(shift(Z_i)^alpha)  [super-linear]
#   CYFER (always exponential) is fit to all four datasets.
#
# METRICS:
#   - Pearson correlation of estimated vs. true cell fate potential.
#   - Jaccard index between top identified features and true causal features.
#   - Gini coefficient of imputed fate potential distribution.
#
# EXPECTED RESULT:
#   CYFER should recover a monotone-equivalent ranking of fate potentials under
#   all smooth, monotone growth models (A-C), since any monotone transformation
#   of Z preserves feature relevance. The power-law (D) also preserves rank
#   when Z_i > 0. Performance degrades gracefully as models diverge from
#   exponential, but the rank-based metrics (Jaccard) should remain high.
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
  2 * sum((seq_len(n)) * x) / (n * sum(x)) - (n + 1) / n
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

n_cells      <- 1500
n_clones     <- 60
n_features   <- 30
n_causal     <- 5          # features with non-zero beta
sigma_between <- 1.0       # between-clone expression variance
sigma_within  <- 0.5       # within-clone expression variance
n_replicates  <- 20        # Monte Carlo replicates

# True coefficients: first n_causal features matter, rest are zero
true_beta <- c(rep(1.0, n_causal), rep(0.0, n_features - n_causal))
names(true_beta) <- paste0("f", seq_len(n_features))
true_causal_features <- names(true_beta[true_beta != 0])

# Cell-to-clone assignment (balanced)
cells_per_clone <- n_cells %/% n_clones
clone_ids       <- rep(seq_len(n_clones), each = cells_per_clone)[seq_len(n_cells)]
clone_labels    <- paste0("clone:", clone_ids)

# ---- Growth model definitions ------------------------------------------------

# (A) Exponential (CYFER's assumed model)
simulate_exponential <- function(Z) {
  rpois(length(Z), lambda = exp(Z))
}

# (B) Linear: counts proportional to shifted/scaled Z (no log link)
simulate_linear <- function(Z) {
  # Shift so minimum lambda = 0.5
  Z_shifted <- Z - min(Z) + 0.5
  rpois(length(Z_shifted), lambda = Z_shifted)
}

# (C) Logistic/saturating: growth bounded by cap = 20 expected progeny
simulate_logistic <- function(Z, cap = 20) {
  lambda <- cap / (1 + exp(-Z))
  rpois(length(lambda), lambda)
}

# (D) Power-law: count ~ Poisson(shift(Z)^alpha), alpha=2
simulate_power <- function(Z, alpha = 2) {
  Z_pos <- Z - min(Z) + 1   # ensure all values > 1 for power law
  lambda <- Z_pos^alpha
  rpois(length(lambda), lambda)
}

growth_models <- list(
  exponential = simulate_exponential,
  linear      = simulate_linear,
  logistic    = simulate_logistic,
  power_law   = simulate_power
)

# ---- Single-run helper -------------------------------------------------------

run_one_replicate <- function(seed_val) {
  set.seed(seed_val)

  # Generate hierarchical feature matrix
  clone_centers <- mvrnorm(n_clones,
                           mu    = rep(0, n_features),
                           Sigma = sigma_between^2 * diag(n_features))
  rownames(clone_centers) <- paste0("clone:", seq_len(n_clones))
  colnames(clone_centers) <- paste0("f", seq_len(n_features))

  X <- t(sapply(seq_len(n_cells), function(i) {
    k <- clone_ids[i]
    clone_centers[k, ] +
      mvrnorm(1, mu = rep(0, n_features),
              Sigma = sigma_within^2 * diag(n_features))
  }))
  rownames(X) <- paste0("cell:", seq_len(n_cells))
  colnames(X) <- paste0("f", seq_len(n_features))
  X <- scale(X)

  # True log fate potential (centered)
  true_Z <- as.numeric(X %*% true_beta)
  names(true_Z) <- rownames(X)
  true_Z <- true_Z - mean(true_Z)   # center at 0

  # Fit CYFER under each growth model
  sapply(names(growth_models), function(model_name) {
    sim_fn <- growth_models[[model_name]]

    # Generate per-cell future counts under this growth model
    cell_future <- sim_fn(true_Z)
    names(cell_future) <- names(true_Z)

    # Aggregate to clone-level (what CYFER actually observes)
    lineage_future_count <- tapply(cell_future, clone_labels, sum)

    # Keep only clones with future count > 0 (others are extinct clones;
    # they are still informative but require count >= 1 for the log term)
    valid_clones <- names(lineage_future_count[lineage_future_count > 0])
    if (length(valid_clones) < 5) return(NA)

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

      # Estimated fate potential on the linear predictor scale
      Z_hat <- as.numeric(X_sub %*% final_fit$coefficient_vec[-1]) +
               final_fit$coefficient_vec[1]
      names(Z_hat) <- rownames(X_sub)

      cor(Z_sub, Z_hat)
    }, error = function(e) NA)
  })
}

# ---- Run all replicates ------------------------------------------------------

message("Running ", n_replicates, " replicates across ", length(growth_models), " growth models...")

replicate_cors <- lapply(seq_len(n_replicates), run_one_replicate)
cor_matrix <- do.call(rbind, replicate_cors)
colnames(cor_matrix) <- names(growth_models)

# ---- Also run single detailed example for Jaccard and Gini ------------------

set.seed(999)

clone_centers_main <- mvrnorm(n_clones, mu = rep(0, n_features),
                              Sigma = sigma_between^2 * diag(n_features))
rownames(clone_centers_main) <- paste0("clone:", seq_len(n_clones))
colnames(clone_centers_main) <- paste0("f", seq_len(n_features))

X_main <- t(sapply(seq_len(n_cells), function(i) {
  k <- clone_ids[i]
  clone_centers_main[k, ] +
    mvrnorm(1, mu = rep(0, n_features),
            Sigma = sigma_within^2 * diag(n_features))
}))
rownames(X_main) <- paste0("cell:", seq_len(n_cells))
colnames(X_main) <- paste0("f", seq_len(n_features))
X_main <- scale(X_main)
true_Z_main <- as.numeric(X_main %*% true_beta)
names(true_Z_main) <- rownames(X_main)
true_Z_main <- true_Z_main - mean(true_Z_main)

detailed_results <- lapply(names(growth_models), function(model_name) {
  sim_fn <- growth_models[[model_name]]
  cell_future <- sim_fn(true_Z_main)
  names(cell_future) <- names(true_Z_main)
  lineage_future_count <- tapply(cell_future, clone_labels, sum)
  valid_clones <- names(lineage_future_count[lineage_future_count > 0])
  keep_idx  <- which(clone_labels %in% valid_clones)
  X_sub     <- X_main[keep_idx, , drop = FALSE]
  Z_sub     <- true_Z_main[keep_idx]
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

    top_k     <- n_causal * 3   # allow generous set
    est_feat  <- top_features_by_cor(X_sub, Z_hat, k = top_k)
    jac_val   <- jaccard(true_causal_features, est_feat[seq_len(n_causal)])
    gini_val  <- gini_coef(exp(Z_hat))
    r_val     <- cor(Z_sub, Z_hat)

    data.frame(model = model_name, correlation = r_val,
               jaccard = jac_val, gini = gini_val,
               n_valid_clones = length(valid_clones),
               stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(model = model_name, correlation = NA, jaccard = NA,
               gini = NA, n_valid_clones = NA, stringsAsFactors = FALSE)
  })
})

detailed_df <- do.call(rbind, detailed_results)

# ---- Print results -----------------------------------------------------------

cat("\n=== Detailed results (single run) ===\n")
print(detailed_df)

cat("\n=== Correlation summary across", n_replicates, "replicates ===\n")
summary_df <- data.frame(
  model  = colnames(cor_matrix),
  mean_r = round(colMeans(cor_matrix, na.rm = TRUE), 3),
  sd_r   = round(apply(cor_matrix, 2, sd, na.rm = TRUE), 3),
  median_r = round(apply(cor_matrix, 2, median, na.rm = TRUE), 3)
)
print(summary_df)

# ---- Save results ------------------------------------------------------------

saveRDS(
  list(
    detailed      = detailed_df,
    replicate_cors = cor_matrix,
    summary        = summary_df,
    params = list(n_cells = n_cells, n_clones = n_clones,
                  n_features = n_features, n_causal = n_causal)
  ),
  file = "sim1_growth_modes_results.rds"
)

message("sim1 complete. Results saved to sim1_growth_modes_results.rds")
