# Writeup21 demo: can the per-clone squared AED span 0 to 2?
# Kevin Z. Lin (drafted by Claude), 2026-09-22
#
# Generation only, no fitting. Two questions. First: with the within-clone
# spread sigma_w varying across clones, does the per-clone squared AED
# (paper's clonal variability score) reach a 0-to-2 range, and what does the
# mean over clones do at the same time? Also records the t2-only Gini each
# setting produces, because clone-varying sigma_w on the causal coordinate
# raises the expected clone size of the spread-out clones (Jensen's
# inequality on the log-normal mean: E[Y_l] has a factor
# exp(beta^2 sigma_l^2 / 2)) and so moves the Gini on its own. Second: if the
# clone-varying part is put on the non-causal coordinates only
# (`spread_variation = "noncausal"`), so that AED and clone size are
# independent by construction, is the AED range preserved and does the
# coupling vanish?
#
# The generator follows simulation_design_claude.md Sections 2.1-2.3 for t1
# cells: hierarchical latent state, one causal coordinate, negative-binomial
# counts on 2,000 genes, PCA to 10 dimensions on log-normalized counts.
# Runs in about four minutes. Writes one CSV and prints the same table.

library(MASS)

rm(list = ls())

# Paths ------------------------------------------------------------------------

repo_dir <- file.path("/Users/kevinlin/Library/CloudStorage/Dropbox",
                      "Collaboration-and-People/archive/Nancy/multiomeFate",
                      "git/multiomeFate_analysis")
csv_dir <- file.path(repo_dir, "csv", "kevin", "Writeup21")
csv_file <- file.path(csv_dir, "demo_aed_range_claude.csv")

# Helpers ----------------------------------------------------------------------

# Gini of a non-negative vector, identical to dineq::gini.wtd on such input.
.gini_coef <- function(x){
  x <- pmax(x, 0)
  x <- sort(x)
  n <- length(x)
  if(n == 0 || sum(x) == 0) return(0)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

# Mean pairwise squared Euclidean distance among the rows of a matrix, without
# forming the pairs: mean_{i != j} ||x_i - x_j||^2 = 2 * sum_k var_k, where
# var_k is the (n - 1)-denominator sample variance of column k.
.mean_sq_dist <- function(mat){
  if(nrow(mat) < 2) return(NA_real_)
  2 * sum(apply(mat, 2, stats::var))
}

# Mean pairwise (un-squared) Euclidean distance, for reference only.
.mean_dist <- function(mat){
  if(nrow(mat) < 2) return(NA_real_)
  mean(stats::dist(mat))
}

# Generate the t1 cells of one dataset and their t2 clone sizes.
#
# The latent total variance per coordinate is tau^2 + mean(sigma_l^2) =
# latent_scale^2, and h2 is the between-clone share of it. AED depends on h2
# only; the Gini depends on latent_scale as well (it sets the spread of Z).
# sigma_l^2 = sigma_w^2 * g_l with g_l ~ Gamma(shape = kappa, rate = kappa),
# mean 1; kappa = Inf gives every clone the same spread.
# spread_variation = "isotropic" applies g_l to all d coordinates;
# "noncausal" applies it to coordinates 2..d only and leaves the causal
# (first) coordinate at the shared sigma_w, so a clone's spread on the
# expansion axis, and hence its expected t2 size, does not depend on g_l.
.generate_t1 <- function(h2,
                         kappa,
                         beta = 1.5,
                         d = 10,
                         latent_scale = 1,
                         num_causal_genes = 100,
                         num_cells_per_clone = 10,
                         num_clones = 100,
                         num_genes = 2000,
                         spread_variation = "isotropic",
                         t2_total = 3000,
                         theta = 10,
                         seed_number = 10){
  stopifnot(h2 >= 0, h2 <= 1, kappa > 0, latent_scale > 0,
            spread_variation %in% c("isotropic", "noncausal"))
  tau <- latent_scale * sqrt(h2)
  sigma_w <- latent_scale * sqrt(1 - h2)

  if(!is.null(seed_number)) set.seed(seed_number)
  num_cells <- num_clones * num_cells_per_clone
  clone_vec <- rep(paste0("clone", seq_len(num_clones)),
                   each = num_cells_per_clone)

  # Clone centres and clone-specific within-clone spreads.
  centre_mat <- MASS::mvrnorm(num_clones, rep(0, d), tau^2 * diag(d))
  if(is.finite(kappa)){
    spread_multiplier_vec <- stats::rgamma(num_clones, shape = kappa,
                                           rate = kappa)
  } else {
    spread_multiplier_vec <- rep(1, num_clones)
  }
  # Per-clone, per-coordinate spread: rows are clones, columns coordinates.
  sigma_mat <- matrix(sigma_w * sqrt(spread_multiplier_vec),
                      nrow = num_clones, ncol = d)
  if(spread_variation == "noncausal") sigma_mat[, 1] <- sigma_w

  clone_idx_vec <- rep(seq_len(num_clones), each = num_cells_per_clone)
  latent_mat <- centre_mat[clone_idx_vec, , drop = FALSE] +
    matrix(stats::rnorm(num_cells * d), nrow = num_cells, ncol = d) *
    sigma_mat[clone_idx_vec, , drop = FALSE]

  # Fate potential on the first coordinate; beta_0 solved so the expected
  # number of t2 cells given these latent states equals t2_total.
  beta_0 <- log(t2_total / sum(exp(beta * latent_mat[, 1])))
  z_vec <- beta_0 + beta * latent_mat[, 1]
  progeny_vec <- stats::rpois(num_cells, exp(z_vec))
  t2_size_vec <- tapply(progeny_vec, clone_vec, sum)[unique(clone_vec)]

  # Gene loadings: a sparse causal column, dense small loadings elsewhere.
  loading_mat <- matrix(stats::rnorm(num_genes * d, sd = 0.25),
                        nrow = num_genes, ncol = d)
  loading_mat[, 1] <- 0
  causal_idx_vec <- seq_len(num_causal_genes)
  loading_mat[causal_idx_vec, 1] <- stats::runif(num_causal_genes, 0.5, 1) *
    rep(c(1, -1), length.out = num_causal_genes)

  # Negative-binomial counts with log-normal baselines and library sizes.
  baseline_vec <- stats::rnorm(num_genes, mean = 0, sd = 1)
  log_mu_mat <- sweep(latent_mat %*% t(loading_mat), 2, baseline_vec, "+")
  mu_mat <- exp(log_mu_mat)
  mu_mat <- mu_mat / rowSums(mu_mat)
  library_vec <- stats::rlnorm(num_cells, meanlog = log(5000), sdlog = 0.3)
  mu_mat <- mu_mat * library_vec
  count_mat <- matrix(stats::rnbinom(length(mu_mat), mu = as.numeric(mu_mat),
                                     size = theta),
                      nrow = num_cells, ncol = num_genes)

  list(clone_vec = clone_vec,
       count_mat = count_mat,
       latent_mat = latent_mat,
       sigma_mat = sigma_mat,
       t2_size_vec = t2_size_vec)
}

# Log-normalize, PCA to 10 dimensions, per-clone squared and un-squared AED.
.compute_aed <- function(count_mat, clone_vec, d = 10){
  lognorm_mat <- log1p(count_mat / rowSums(count_mat) * 1e4)
  keep_gene_vec <- which(apply(lognorm_mat, 2, stats::var) > 0)
  pca_res <- stats::prcomp(lognorm_mat[, keep_gene_vec, drop = FALSE],
                           rank. = d, center = TRUE, scale. = FALSE)
  embedding_mat <- pca_res$x[, seq_len(d), drop = FALSE]

  sq_random <- .mean_sq_dist(embedding_mat)
  dist_random <- .mean_dist(embedding_mat)
  clone_name_vec <- unique(clone_vec)
  aed_sq_vec <- sapply(clone_name_vec, function(clone){
    .mean_sq_dist(embedding_mat[clone_vec == clone, , drop = FALSE]) / sq_random
  })
  aed_unsq_vec <- sapply(clone_name_vec, function(clone){
    .mean_dist(embedding_mat[clone_vec == clone, , drop = FALSE]) / dist_random
  })

  list(aed_sq_vec = aed_sq_vec,
       aed_unsq_vec = aed_unsq_vec)
}

# The grid ---------------------------------------------------------------------

# h2 = tau^2 / (tau^2 + mean sigma_l^2): expected squared AED is about 1 - h2.
# kappa: Inf is a shared sigma_w; smaller is more clone-to-clone variation.
grid_a_df <- expand.grid(h2 = c(0.9, 0.75, 0.5, 0.25, 0.1, 0),
                         kappa = c(Inf, 3, 1.5, 1),
                         latent_scale = 1,
                         seed_number = c(10, 20),
                         spread_variation = "isotropic")
# Second grid: at the top AED level, does shrinking the latent scale bring the
# Gini down (and does the count-noise floor of the AED stay small)?
grid_b_df <- expand.grid(h2 = 0,
                         kappa = c(Inf, 1.5),
                         latent_scale = c(0.3, 0.5, 0.7),
                         seed_number = c(10, 20),
                         spread_variation = "isotropic")
# Third grid: the same levels with the clone-varying spread on the non-causal
# coordinates only. kappa = Inf is identical in both modes, so it is not
# repeated.
grid_c_df <- expand.grid(h2 = c(0.9, 0.75, 0.5, 0.25, 0.1, 0),
                         kappa = c(3, 1.5, 1),
                         latent_scale = 1,
                         seed_number = c(10, 20),
                         spread_variation = "noncausal")
grid_d_df <- expand.grid(h2 = 0,
                         kappa = 1.5,
                         latent_scale = c(0.3, 0.5, 0.7),
                         seed_number = c(10, 20),
                         spread_variation = "noncausal")
grid_df <- rbind(grid_a_df, grid_b_df, grid_c_df, grid_d_df)
grid_df$spread_variation <- as.character(grid_df$spread_variation)

result_list <- vector("list", nrow(grid_df))
for(i in seq_len(nrow(grid_df))){
  h2 <- grid_df$h2[i]
  kappa <- grid_df$kappa[i]
  latent_scale <- grid_df$latent_scale[i]
  seed_number <- grid_df$seed_number[i]
  spread_variation <- grid_df$spread_variation[i]
  print(paste0(Sys.time(), " | setting ", i, " of ", nrow(grid_df),
               ": h2 = ", h2, ", kappa = ", kappa, ", scale = ", latent_scale,
               ", seed = ", seed_number, ", spread = ", spread_variation))

  data_list <- .generate_t1(h2 = h2,
                            kappa = kappa,
                            latent_scale = latent_scale,
                            spread_variation = spread_variation,
                            seed_number = seed_number)
  aed_list <- .compute_aed(data_list$count_mat, data_list$clone_vec)
  aed_sq_vec <- aed_list$aed_sq_vec

  result_list[[i]] <- data.frame(
    h2 = h2,
    kappa = kappa,
    latent_scale = latent_scale,
    seed = seed_number,
    spread_variation = spread_variation,
    aed_sq_mean = mean(aed_sq_vec),
    aed_sq_min = min(aed_sq_vec),
    aed_sq_q05 = unname(stats::quantile(aed_sq_vec, 0.05)),
    aed_sq_median = stats::median(aed_sq_vec),
    aed_sq_q95 = unname(stats::quantile(aed_sq_vec, 0.95)),
    aed_sq_max = max(aed_sq_vec),
    aed_unsq_mean = mean(aed_list$aed_unsq_vec),
    gini_t2 = .gini_coef(data_list$t2_size_vec),
    t2_total = sum(data_list$t2_size_vec),
    t2_max_clone = max(data_list$t2_size_vec),
    t2_zero_clones = sum(data_list$t2_size_vec == 0),
    spearman_aed_t2size = stats::cor(aed_sq_vec, data_list$t2_size_vec,
                                     method = "spearman"))
}
result_df <- do.call(rbind, result_list)

# Output -----------------------------------------------------------------------

dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(result_df, csv_file, row.names = FALSE)

numeric_col_vec <- sapply(result_df, is.numeric)
result_df[numeric_col_vec] <- round(result_df[numeric_col_vec], 2)
print(result_df)
print(paste0("Written to ", csv_file))
