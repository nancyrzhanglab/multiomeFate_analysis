# func_generate_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-22
#
# The generator behind both Figure 4 sweeps (simulation_design_claude.md,
# Sections 2 and 5). Three layers:
#   .generate_latent()   t1 latent states, fate potential, t2 clone sizes
#                        (milliseconds; what the Gini calibration loops over)
#   generate_dataset()   the full dataset: t1 cells plus their t2 children,
#                        gene loadings, negative-binomial counts for every cell
#   calibrate_scale()    bisects the latent scale `s` so the t2 Gini hits a
#                        target at a given h2 (Section 5)
# plus the two statistics the axes are defined by, gini_coef() and
# mean_sq_dist().
#
# Conventions. `h2` is the between-clone share of latent variance, so the
# clone-centre spread is tau = s * sqrt(h2) and the typical within-clone
# spread is sigma_w = s * sqrt(1 - h2). The clone-varying multiplier
# g_l ~ Gamma(kappa, kappa) (mean 1) scales the within-clone variance of
# either every coordinate (`spread_variation = "isotropic"`) or of the
# non-causal coordinates 2..d only (`"noncausal"`, the sweeps' choice). The
# first latent coordinate is the expansion axis: Z_i = beta_0 + beta * s_i1
# and each t1 cell leaves Poisson(exp(Z_i)) children at t2.

# Statistics ------------------------------------------------------------------

#' Gini coefficient of a non-negative vector
#'
#' Identical to `dineq::gini.wtd` on non-negative input, which is what the
#' paper's real-data numbers use.
#'
#' @param x numeric vector; negatives are clipped to 0.
#'
#' @returns numeric scalar in [0, 1 - 1/length(x)]; 0 when the vector sums to 0.
gini_coef <- function(x){
  x <- pmax(x, 0)
  x <- sort(x)
  n <- length(x)
  if(n == 0 || sum(x) == 0) return(0)
  2 * sum(seq_len(n) * x) / (n * sum(x)) - (n + 1) / n
}

#' Mean pairwise squared Euclidean distance among the rows of a matrix
#'
#' Computed without forming the pairs: mean_{i != j} ||x_i - x_j||^2 equals
#' twice the sum of the per-column sample variances (n - 1 denominator).
#'
#' @param mat numeric matrix, rows are points.
#'
#' @returns numeric scalar; `NA` when there are fewer than two rows.
mean_sq_dist <- function(mat){
  if(nrow(mat) < 2) return(NA_real_)
  2 * sum(apply(mat, 2, stats::var))
}

# Latent layer ----------------------------------------------------------------

#' Draw the t1 latent states, fate potentials and t2 clone sizes
#'
#' The cheap core of the generator, shared by `generate_dataset()` and the
#' calibration. No genes, no counts.
#'
#' @param h2 between-clone share of the latent variance, in [0, 1].
#' @param latent_scale the overall latent scale `s`, positive.
#' @param beta effect of the causal coordinate on the log fate potential.
#' @param d latent dimension.
#' @param kappa shape of the Gamma(kappa, kappa) clone-varying spread
#'   multiplier; `Inf` gives every clone the same spread.
#' @param num_cells_per_clone t1 cells per clone (equal across clones).
#' @param num_clones number of clones.
#' @param spread_variation `"noncausal"` or `"isotropic"`; where the
#'   clone-varying multiplier acts (see file header).
#' @param t2_total expected total number of t2 cells; `beta_0` is solved for it.
#' @param seed_number seed, `NULL` to leave the stream alone.
#'
#' @returns a list with `beta_0`, `centre_mat` (clones by d), `clone_vec`
#'   (character, one per t1 cell), `latent_mat` (t1 cells by d),
#'   `progeny_vec` (children per t1 cell), `sigma_mat` (clones by d, the
#'   per-coordinate within-clone standard deviations),
#'   `spread_multiplier_vec` (g_l per clone), `t2_size_vec` (named by clone,
#'   zeros kept) and `z_vec` (the true log fate potential per t1 cell).
#' @noRd
.generate_latent <- function(h2,
                             latent_scale,
                             beta = 1.5,
                             d = 10,
                             kappa = 1.5,
                             num_cells_per_clone = 10,
                             num_clones = 100,
                             spread_variation = "noncausal",
                             t2_total = 3000,
                             seed_number = 10){
  tau <- latent_scale * sqrt(h2)
  sigma_w <- latent_scale * sqrt(1 - h2)
  if(!is.null(seed_number)) set.seed(seed_number)

  num_cells <- num_clones * num_cells_per_clone
  clone_name_vec <- paste0("clone", seq_len(num_clones))
  clone_idx_vec <- rep(seq_len(num_clones), each = num_cells_per_clone)
  clone_vec <- clone_name_vec[clone_idx_vec]

  centre_mat <- matrix(stats::rnorm(num_clones * d, sd = tau),
                       nrow = num_clones, ncol = d)
  if(is.finite(kappa)){
    spread_multiplier_vec <- stats::rgamma(num_clones, shape = kappa,
                                           rate = kappa)
  } else {
    spread_multiplier_vec <- rep(1, num_clones)
  }
  sigma_mat <- matrix(sigma_w * sqrt(spread_multiplier_vec),
                      nrow = num_clones, ncol = d)
  if(spread_variation == "noncausal") sigma_mat[, 1] <- sigma_w

  latent_mat <- centre_mat[clone_idx_vec, , drop = FALSE] +
    matrix(stats::rnorm(num_cells * d), nrow = num_cells, ncol = d) *
    sigma_mat[clone_idx_vec, , drop = FALSE]

  # beta_0 makes the expected number of t2 cells, given these latent states,
  # equal t2_total at every level of both sweeps.
  beta_0 <- log(t2_total / sum(exp(beta * latent_mat[, 1])))
  z_vec <- beta_0 + beta * latent_mat[, 1]
  progeny_vec <- stats::rpois(num_cells, exp(z_vec))
  t2_size_vec <- as.numeric(tapply(progeny_vec, clone_vec, sum)[clone_name_vec])
  names(t2_size_vec) <- clone_name_vec

  list(beta_0 = beta_0,
       centre_mat = centre_mat,
       clone_vec = clone_vec,
       latent_mat = latent_mat,
       progeny_vec = progeny_vec,
       sigma_mat = sigma_mat,
       spread_multiplier_vec = spread_multiplier_vec,
       t2_size_vec = t2_size_vec,
       z_vec = z_vec)
}

# Full dataset ----------------------------------------------------------------

#' Generate one dataset: t1 cells, their t2 children, and counts for all
#'
#' Sections 2.1 to 2.4 of the design memo. The t2 population is the multiset
#' of children of the t1 cells (no capture subsampling), each child's latent
#' state is a rho-weighted mix of its parent and the clone centre plus fresh
#' within-clone noise, a shared shift `delta` and a clone-specific shift
#' `delta_l`, both on the non-causal coordinates. Counts for every cell come
#' from one gene model.
#'
#' @param h2 between-clone share of the latent variance, in [0, 1].
#' @param latent_scale the overall latent scale `s`, positive.
#' @param beta effect of the causal coordinate on the log fate potential.
#' @param d latent dimension.
#' @param delta_multiplier norm of the shared t2 shift in units of the
#'   latent scale: `||delta|| = delta_multiplier * latent_scale`, so that t1
#'   and t2 cells barely overlap in the top PCs.
#' @param kappa shape of the clone-varying spread multiplier; `Inf` for none.
#' @param library_meanlog,library_sdlog log-normal library-size parameters.
#' @param loading_sd standard deviation of the dense non-causal loadings.
#' @param num_causal_genes genes loading on the causal coordinate, half with
#'   positive and half with negative loading.
#' @param num_cells_per_clone t1 cells per clone.
#' @param num_clones number of clones.
#' @param num_genes number of genes.
#' @param rho child-parent latent correlation at t2.
#' @param spread_variation `"noncausal"` or `"isotropic"`.
#' @param t2_total expected total number of t2 cells.
#' @param tau_delta standard deviation of the clone-specific t2 shift
#'   `delta_l` on each non-causal coordinate.
#' @param theta negative-binomial size (dispersion) parameter.
#' @param seed_number seed, `NULL` to leave the stream alone. The three
#'   stochastic stages use `10 * seed_number + 0, 1, 2`, so consecutive
#'   seeds share no stream and no dataset reuses a calibration draw (which
#'   seed `.generate_latent()` with the plain seed).
#'
#' @returns a list with `cell_df` (one row per cell in the row order of
#'   `count_mat`: `cell_id`, `time_info` (`"t1"`/`"t2"`), `clone_id`,
#'   `parent_id` (`NA` for t1 cells)), `count_mat` (integer matrix, all cells
#'   by genes, with row and column names), `latent_mat` (all cells by d),
#'   `loading_mat` (genes by d), `param_list` (every argument), `t1_idx`
#'   and `t2_idx` (row indices), `t2_size_vec` (named by clone, zeros kept),
#'   `z_true_vec` (named log fate potential of the t1 cells) and `latent`
#'   (the full `.generate_latent()` output).
generate_dataset <- function(h2,
                             latent_scale,
                             beta = 1.5,
                             d = 10,
                             delta_multiplier = 3,
                             kappa = 1.5,
                             library_meanlog = log(5000),
                             library_sdlog = 0.3,
                             loading_sd = 0.25,
                             num_causal_genes = 100,
                             num_cells_per_clone = 10,
                             num_clones = 100,
                             num_genes = 2000,
                             rho = 0.8,
                             spread_variation = "noncausal",
                             t2_total = 3000,
                             tau_delta = 0.6,
                             theta = 10,
                             seed_number = 10){
  stopifnot(h2 >= 0, h2 <= 1, latent_scale > 0, kappa > 0, d >= 2,
            rho >= 0, rho <= 1, tau_delta >= 0, delta_multiplier >= 0,
            num_causal_genes <= num_genes,
            spread_variation %in% c("isotropic", "noncausal"))
  param_list <- list(beta = beta, d = d, delta_multiplier = delta_multiplier,
                     h2 = h2, kappa = kappa, latent_scale = latent_scale,
                     library_meanlog = library_meanlog,
                     library_sdlog = library_sdlog, loading_sd = loading_sd,
                     num_causal_genes = num_causal_genes,
                     num_cells_per_clone = num_cells_per_clone,
                     num_clones = num_clones, num_genes = num_genes,
                     rho = rho, seed_number = seed_number,
                     spread_variation = spread_variation, t2_total = t2_total,
                     tau_delta = tau_delta, theta = theta)

  # Step 1: t1 latent states, fate potentials and progeny counts.
  latent <- .generate_latent(h2 = h2,
                             latent_scale = latent_scale,
                             beta = beta,
                             d = d,
                             kappa = kappa,
                             num_cells_per_clone = num_cells_per_clone,
                             num_clones = num_clones,
                             spread_variation = spread_variation,
                             t2_total = t2_total,
                             seed_number = if(is.null(seed_number)) NULL else
                               10 * seed_number)
  num_t1 <- nrow(latent$latent_mat)
  clone_idx_vec <- rep(seq_len(num_clones), each = num_cells_per_clone)

  # Step 2: the t2 population is the multiset of children.
  # `rep()` on the t1 indices expands each parent into its progeny_vec[i]
  # children; a parent with 0 children contributes no row.
  if(!is.null(seed_number)) set.seed(10 * seed_number + 1)
  parent_idx_vec <- rep(seq_len(num_t1), times = latent$progeny_vec)
  num_t2 <- length(parent_idx_vec)
  noncausal_idx_vec <- seq(2, d)

  # Shared shift: a fixed direction on the non-causal coordinates, scaled to
  # delta_multiplier standard deviations of the t1 cloud (per coordinate the
  # t1 latent standard deviation is latent_scale).
  delta_vec <- rep(0, d)
  direction_vec <- stats::rnorm(d - 1)
  delta_vec[noncausal_idx_vec] <- direction_vec / sqrt(sum(direction_vec^2)) *
    delta_multiplier * latent_scale
  # Clone-specific shift on the non-causal coordinates, one draw per clone.
  clone_shift_mat <- matrix(0, nrow = num_clones, ncol = d)
  clone_shift_mat[, noncausal_idx_vec] <- stats::rnorm(num_clones * (d - 1),
                                                       sd = tau_delta)

  child_latent_mat <- matrix(0, nrow = num_t2, ncol = d)
  if(num_t2 > 0){
    child_clone_idx_vec <- clone_idx_vec[parent_idx_vec]
    child_latent_mat <- rho * latent$latent_mat[parent_idx_vec, , drop = FALSE] +
      (1 - rho) * latent$centre_mat[child_clone_idx_vec, , drop = FALSE] +
      sqrt(1 - rho^2) * latent$sigma_mat[child_clone_idx_vec, , drop = FALSE] *
      matrix(stats::rnorm(num_t2 * d), nrow = num_t2, ncol = d) +
      matrix(delta_vec, nrow = num_t2, ncol = d, byrow = TRUE) +
      clone_shift_mat[child_clone_idx_vec, , drop = FALSE]
  }
  latent_mat <- rbind(latent$latent_mat, child_latent_mat)
  num_cells <- nrow(latent_mat)

  # Step 3: gene loadings, a sparse causal column and dense small loadings
  # on the other coordinates, then negative-binomial counts for every cell.
  if(!is.null(seed_number)) set.seed(10 * seed_number + 2)
  loading_mat <- matrix(stats::rnorm(num_genes * d, sd = loading_sd),
                        nrow = num_genes, ncol = d)
  loading_mat[, 1] <- 0
  causal_idx_vec <- seq_len(num_causal_genes)
  loading_mat[causal_idx_vec, 1] <- stats::runif(num_causal_genes, 0.5, 1) *
    rep(c(1, -1), length.out = num_causal_genes)
  baseline_vec <- stats::rnorm(num_genes, mean = 0, sd = 1)

  log_mu_mat <- sweep(latent_mat %*% t(loading_mat), 2, baseline_vec, "+")
  mu_mat <- exp(log_mu_mat)
  mu_mat <- mu_mat / rowSums(mu_mat)
  library_vec <- stats::rlnorm(num_cells, meanlog = library_meanlog,
                               sdlog = library_sdlog)
  mu_mat <- mu_mat * library_vec
  count_mat <- matrix(stats::rnbinom(length(mu_mat), mu = as.numeric(mu_mat),
                                     size = theta),
                      nrow = num_cells, ncol = num_genes)

  cell_id_vec <- c(paste0("t1_cell", seq_len(num_t1)),
                   if(num_t2 > 0) paste0("t2_cell", seq_len(num_t2)) else NULL)
  rownames(count_mat) <- cell_id_vec
  colnames(count_mat) <- paste0("gene", seq_len(num_genes))
  rownames(latent_mat) <- cell_id_vec

  cell_df <- data.frame(
    cell_id = cell_id_vec,
    time_info = c(rep("t1", num_t1), rep("t2", num_t2)),
    clone_id = c(latent$clone_vec, latent$clone_vec[parent_idx_vec]),
    parent_id = c(rep(NA_character_, num_t1),
                  cell_id_vec[parent_idx_vec]),
    stringsAsFactors = FALSE)
  z_true_vec <- latent$z_vec
  names(z_true_vec) <- cell_id_vec[seq_len(num_t1)]

  stopifnot(sum(latent$t2_size_vec) == num_t2,
            nrow(count_mat) == nrow(cell_df))

  list(cell_df = cell_df,
       count_mat = count_mat,
       latent = latent,
       latent_mat = latent_mat,
       loading_mat = loading_mat,
       param_list = param_list,
       t1_idx = seq_len(num_t1),
       t2_idx = if(num_t2 > 0) num_t1 + seq_len(num_t2) else integer(0),
       t2_size_vec = latent$t2_size_vec,
       z_true_vec = z_true_vec)
}

# Calibration -----------------------------------------------------------------

#' Mean realized t2 Gini at a latent scale, over quick latent-only draws
#'
#' @param latent_scale the scale `s` to evaluate.
#' @param h2 between-clone share of the latent variance.
#' @param num_draws number of latent draws to average over.
#' @param seed_number base seed; draw `k` uses `seed_number + k`.
#' @param ... passed to `.generate_latent()`.
#'
#' @returns numeric scalar, the mean t2 Gini.
#' @noRd
.mean_gini_at_scale <- function(latent_scale,
                                h2,
                                num_draws = 10,
                                seed_number = 10,
                                ...){
  gini_vec <- sapply(seq_len(num_draws), function(k){
    latent <- .generate_latent(h2 = h2,
                               latent_scale = latent_scale,
                               seed_number = seed_number + k,
                               ...)
    gini_coef(latent$t2_size_vec)
  })
  mean(gini_vec)
}

#' Bisect the latent scale so the mean t2 Gini hits a target
#'
#' Section 5 of the design memo. The mean squared AED is about `1 - h2`
#' regardless of the scale, and the t2 Gini rises monotonically with the
#' scale at fixed `h2`, so each axis is a one-knob calibration: `h2` is set
#' from the AED target directly and `s` is bisected for the Gini target.
#'
#' @param h2 between-clone share of the latent variance.
#' @param gini_target target mean t2 Gini, above the Poisson floor.
#' @param max_iter maximum bisection steps.
#' @param num_draws latent draws averaged per candidate scale.
#' @param scale_range initial bracket on `s`.
#' @param tol stop when the mean Gini is within `tol` of the target.
#' @param seed_number base seed for the draws.
#' @param verbose numeric.
#' @param ... passed to `.generate_latent()` (kappa, spread_variation, ...).
#'
#' @returns a list with `gini_realized` (mean over the draws at the chosen
#'   scale), `latent_scale`, `num_iter`, `sigma_w` and `tau`.
calibrate_scale <- function(h2,
                            gini_target,
                            max_iter = 40,
                            num_draws = 10,
                            scale_range = c(0.01, 8),
                            tol = 0.01,
                            seed_number = 10,
                            verbose = 0,
                            ...){
  stopifnot(h2 >= 0, h2 <= 1, gini_target > 0, gini_target < 1,
            length(scale_range) == 2, scale_range[1] < scale_range[2])

  lower_val <- scale_range[1]
  upper_val <- scale_range[2]
  gini_lower <- .mean_gini_at_scale(lower_val, h2 = h2, num_draws = num_draws,
                                    seed_number = seed_number, ...)
  gini_upper <- .mean_gini_at_scale(upper_val, h2 = h2, num_draws = num_draws,
                                    seed_number = seed_number, ...)
  if(gini_target < gini_lower || gini_target > gini_upper){
    stop("`gini_target` (", gini_target, ") is outside the Gini range [",
         round(gini_lower, 3), ", ", round(gini_upper, 3),
         "] reachable over `scale_range`")
  }

  num_iter <- 0
  mid_val <- NA_real_
  gini_mid <- NA_real_
  for(iter in seq_len(max_iter)){
    num_iter <- iter
    mid_val <- sqrt(lower_val * upper_val)
    gini_mid <- .mean_gini_at_scale(mid_val, h2 = h2, num_draws = num_draws,
                                    seed_number = seed_number, ...)
    if(verbose > 1){
      print(paste0("  bisection ", iter, ": s = ", signif(mid_val, 4),
                   ", mean Gini = ", round(gini_mid, 4)))
    }
    if(abs(gini_mid - gini_target) <= tol) break
    if(gini_mid < gini_target){
      lower_val <- mid_val
    } else {
      upper_val <- mid_val
    }
  }
  if(verbose > 0){
    print(paste0("calibrate_scale: h2 = ", h2, ", Gini target ", gini_target,
                 " -> s = ", signif(mid_val, 4), " (realized ",
                 round(gini_mid, 3), ", ", num_iter, " iterations)"))
  }

  list(gini_realized = gini_mid,
       latent_scale = mid_val,
       num_iter = num_iter,
       sigma_w = mid_val * sqrt(1 - h2),
       tau = mid_val * sqrt(h2))
}
