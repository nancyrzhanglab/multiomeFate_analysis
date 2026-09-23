# func_methods_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-22
#
# Everything downstream of the generator for the Writeup21 sweeps
# (simulation_design_claude.md, Sections 4.1, 6 and 7): the shared PCA
# embedding, the per-clone squared AED, the operational truth vector, the
# three methods (CYFER, lineage-DE, CoSPAR) each returning a per-gene
# association vector with a p-value, and the metric that scores such a
# vector against the truth.
#
# Every per-gene statistic is signed and scale-free: a Spearman correlation
# with a fate score (CYFER, CoSPAR, truth) or the rank-biserial correlation
# from the Wilcoxon rank-sum test (lineage-DE). The headline metric is the
# Spearman correlation between a method's vector and the truth vector over
# all genes.
#
# Requires: `cospar_flat_io.R` (sourced by the driver) for the CoSPAR
# export/import, and the `multiomeFate` package.

# Embedding and AED -----------------------------------------------------------

#' Log-normalize counts and embed all cells in a PCA fitted on the t1 cells
#'
#' The one embedding every method and the AED share (Section 6.2). PCA is
#' fitted on the t1 cells; t2 cells are projected with the t1 centre and
#' rotation. Genes with zero variance across the t1 cells are dropped from
#' the PCA.
#'
#' @param count_mat integer matrix, all cells by genes, with row names.
#' @param t1_idx integer vector of the t1 rows.
#' @param d number of PCs.
#' @param scale_factor library-size scale before `log1p`.
#'
#' @returns a list with `keep_gene_vec` (genes used in the PCA),
#'   `lognorm_mat` (all cells by genes), `pca_mat` (all cells by d, row
#'   names kept) and `sdev_vec` (the t1 PCA standard deviations).
compute_embedding <- function(count_mat,
                              t1_idx,
                              d = 10,
                              scale_factor = 1e4){
  stopifnot(is.matrix(count_mat), !is.null(rownames(count_mat)),
            length(t1_idx) > d)
  lognorm_mat <- log1p(count_mat / rowSums(count_mat) * scale_factor)
  var_vec <- apply(lognorm_mat[t1_idx, , drop = FALSE], 2, stats::var)
  keep_gene_vec <- which(var_vec > 0)

  pca_res <- stats::prcomp(lognorm_mat[t1_idx, keep_gene_vec, drop = FALSE],
                           rank. = d, center = TRUE, scale. = FALSE)
  pca_mat <- sweep(lognorm_mat[, keep_gene_vec, drop = FALSE], 2,
                   pca_res$center, "-") %*% pca_res$rotation[, seq_len(d)]
  rownames(pca_mat) <- rownames(count_mat)
  colnames(pca_mat) <- paste0("PC", seq_len(d))

  list(keep_gene_vec = keep_gene_vec,
       lognorm_mat = lognorm_mat,
       pca_mat = pca_mat,
       sdev_vec = pca_res$sdev[seq_len(d)])
}

#' Per-clone squared AED in an embedding
#'
#' The paper's clonal variability score, squared (Section 4.1): the mean
#' pairwise squared distance among a clone's cells over the mean pairwise
#' squared distance among all cells, both in the same embedding.
#'
#' @param embedding_mat numeric matrix, t1 cells by dimensions.
#' @param clone_vec character vector of clone labels, one per row.
#'
#' @returns named numeric vector, one AED per clone (`NA` for a clone with
#'   fewer than two cells).
compute_aed <- function(embedding_mat, clone_vec){
  stopifnot(nrow(embedding_mat) == length(clone_vec))
  sq_random <- mean_sq_dist(embedding_mat)
  clone_name_vec <- unique(clone_vec)
  aed_vec <- sapply(clone_name_vec, function(clone){
    mean_sq_dist(embedding_mat[clone_vec == clone, , drop = FALSE]) / sq_random
  })
  names(aed_vec) <- clone_name_vec
  aed_vec
}

#' Between-clone share of the variance of an embedding
#'
#' The `.anova_percentage`-style heritability of the PCA embedding
#' (Section 7): per PC, the between-clone sum of squares over the total,
#' then a variance-weighted mean over PCs.
#'
#' @param embedding_mat numeric matrix, t1 cells by dimensions.
#' @param clone_vec character vector of clone labels, one per row.
#'
#' @returns numeric scalar in [0, 1].
compute_pca_heritability <- function(embedding_mat, clone_vec){
  stopifnot(nrow(embedding_mat) == length(clone_vec))
  clone_factor <- factor(clone_vec)
  total_ss_vec <- apply(embedding_mat, 2, function(x) sum((x - mean(x))^2))
  between_ss_vec <- apply(embedding_mat, 2, function(x){
    clone_mean_vec <- tapply(x, clone_factor, mean)
    clone_size_vec <- tabulate(clone_factor)
    sum(clone_size_vec * (clone_mean_vec - mean(x))^2)
  })
  sum(between_ss_vec) / sum(total_ss_vec)
}

# Per-gene statistics ---------------------------------------------------------

#' Spearman correlation of every gene with a per-cell score
#'
#' Vectorized: rank every column and the score, then one `stats::cor()`
#' call. The p-value is the t-approximation `t = r sqrt((n - 2) / (1 - r^2))`
#' on n - 2 degrees of freedom, which is what `stats::cor.test(method =
#' "spearman")` falls back to in the presence of ties (counts always have
#' ties), so the two agree.
#'
#' @param expr_mat numeric matrix, cells by genes.
#' @param score_vec numeric vector, one per cell.
#'
#' @returns a data frame with `gene`, `stat` (the Spearman correlation) and
#'   `pvalue`, one row per gene in column order.
gene_spearman <- function(expr_mat, score_vec){
  stopifnot(nrow(expr_mat) == length(score_vec))
  n <- nrow(expr_mat)
  rank_mat <- apply(expr_mat, 2, rank)
  rank_score_vec <- rank(score_vec)
  stat_vec <- as.numeric(stats::cor(rank_mat, rank_score_vec))
  # A constant gene has an undefined correlation; call it 0 with p = 1.
  stat_vec[is.na(stat_vec)] <- 0
  t_vec <- stat_vec * sqrt((n - 2) / pmax(1 - stat_vec^2, 1e-12))
  pvalue_vec <- 2 * stats::pt(abs(t_vec), df = n - 2, lower.tail = FALSE)

  data.frame(gene = colnames(expr_mat),
             stat = stat_vec,
             pvalue = pvalue_vec,
             stringsAsFactors = FALSE)
}

#' Rank-biserial correlation of every gene between two cell groups
#'
#' The lineage-DE statistic (Section 6.3): `m_g = 2 U_g / (n_H n_L) - 1`
#' where `U_g` is the Mann-Whitney U of the high group, computed for all
#' genes at once from average ranks, `U_g = sum_{i in H} rank_g[i] -
#' n_H (n_H + 1) / 2`. `U_g` is exactly the `W` that
#' `stats::wilcox.test(x_H, x_L)$statistic` reports. The p-value is the
#' normal approximation with the tie correction and continuity correction
#' that `stats::wilcox.test()` uses when either group exceeds 49 cells.
#'
#' @param expr_mat numeric matrix, cells by genes.
#' @param bool_high_vec logical vector, one per cell, `TRUE` for the high group.
#'
#' @returns a data frame with `gene`, `stat` (the rank-biserial correlation,
#'   in [-1, 1]), `u` (the U statistic) and `pvalue`.
gene_rank_biserial <- function(expr_mat, bool_high_vec){
  stopifnot(nrow(expr_mat) == length(bool_high_vec), is.logical(bool_high_vec))
  num_high <- sum(bool_high_vec)
  num_low <- sum(!bool_high_vec)
  stopifnot(num_high > 0, num_low > 0)
  n <- num_high + num_low

  rank_mat <- apply(expr_mat, 2, rank)
  u_vec <- colSums(rank_mat[bool_high_vec, , drop = FALSE]) -
    num_high * (num_high + 1) / 2
  stat_vec <- 2 * u_vec / (num_high * num_low) - 1

  # Normal approximation of stats::wilcox.test(): variance with the tie
  # term sum(t^3 - t) over tie groups, continuity correction of 1/2.
  tie_term_vec <- apply(rank_mat, 2, function(r){
    tie_vec <- table(r)
    sum(tie_vec^3 - tie_vec)
  })
  var_vec <- num_high * num_low / 12 *
    ((n + 1) - tie_term_vec / (n * (n - 1)))
  z_vec <- u_vec - num_high * num_low / 2
  z_vec <- (z_vec - sign(z_vec) * 0.5) / sqrt(pmax(var_vec, 1e-12))
  pvalue_vec <- 2 * stats::pnorm(abs(z_vec), lower.tail = FALSE)
  pvalue_vec[var_vec <= 0] <- 1

  data.frame(gene = colnames(expr_mat),
             stat = stat_vec,
             u = u_vec,
             pvalue = pvalue_vec,
             stringsAsFactors = FALSE)
}

# Truth and metric ------------------------------------------------------------

#' The operational truth: every gene's Spearman correlation with the true
#' fate potential over the t1 cells
#'
#' @param lognorm_t1_mat log-normalized expression, t1 cells by genes.
#' @param z_true_vec true log fate potential, one per t1 cell.
#' @param q_threshold BH threshold defining the truth *set* for Jaccard.
#'
#' @returns a list with `bool_set_vec` (logical per gene, the truth set) and
#'   `truth_df` (`gene`, `stat`, `pvalue`, `qvalue`).
compute_truth <- function(lognorm_t1_mat,
                          z_true_vec,
                          q_threshold = 0.05){
  truth_df <- gene_spearman(lognorm_t1_mat, z_true_vec)
  truth_df$qvalue <- stats::p.adjust(truth_df$pvalue, method = "BH")

  list(bool_set_vec = truth_df$qvalue < q_threshold,
       truth_df = truth_df)
}

#' Score a method's per-gene vector against the truth
#'
#' @param stat_vec the method's per-gene statistic, in the truth's gene order.
#' @param pvalue_vec the method's per-gene p-value, same order.
#' @param truth_list output of `compute_truth()`.
#' @param q_threshold BH threshold for the method's called set.
#'
#' @returns a one-row data frame with `spearman` (the headline metric),
#'   `pearson`, `jaccard` (called set against the truth set at
#'   `q_threshold`), `num_called` and `num_truth`.
score_gene_stats <- function(stat_vec,
                             pvalue_vec,
                             truth_list,
                             q_threshold = 0.05){
  truth_vec <- truth_list$truth_df$stat
  stopifnot(length(stat_vec) == length(truth_vec),
            length(pvalue_vec) == length(truth_vec))
  if(all(is.na(stat_vec))){
    return(data.frame(spearman = NA_real_, pearson = NA_real_,
                      jaccard = NA_real_, num_called = NA_integer_,
                      num_truth = sum(truth_list$bool_set_vec)))
  }
  bool_called_vec <- stats::p.adjust(pvalue_vec, method = "BH") < q_threshold
  bool_truth_vec <- truth_list$bool_set_vec
  union_size <- sum(bool_called_vec | bool_truth_vec)
  jaccard <- if(union_size == 0) NA_real_ else {
    sum(bool_called_vec & bool_truth_vec) / union_size
  }

  data.frame(spearman = stats::cor(stat_vec, truth_vec, method = "spearman"),
             pearson = stats::cor(stat_vec, truth_vec),
             jaccard = jaccard,
             num_called = sum(bool_called_vec),
             num_truth = sum(bool_truth_vec))
}

# Methods ---------------------------------------------------------------------

#' High clones: the lineage-DE split, shared with CoSPAR's fate labels
#'
#' A clone is "high" when its t2 size exceeds the mean t2 size (Section
#' 6.3): the mean rather than the median, so that a few dominant clones
#' leave a small high group.
#'
#' @param t2_size_vec named numeric vector of t2 clone sizes, zeros included.
#'
#' @returns character vector of the high clone names.
high_clones <- function(t2_size_vec){
  names(t2_size_vec)[t2_size_vec > mean(t2_size_vec)]
}

#' CYFER: fit on the t1 PCs with the zero-count clones kept, then correlate
#' every gene with the estimated fate potential
#'
#' Fails soft: a fitting error returns `NA` statistics and
#' `bool_converged = FALSE` rather than stopping the sweep.
#'
#' @param pca_t1_mat numeric matrix, t1 cells by PCs, with row names.
#' @param clone_vec character vector, the clone of every t1 cell.
#' @param t2_size_vec named numeric vector of t2 clone sizes, zeros kept.
#' @param lognorm_t1_mat log-normalized expression, t1 cells by genes.
#' @param lambda_initial,lambda_sequence_length,num_folds passed to
#'   `multiomeFate::cyfer()`.
#' @param seed_number passed to `cyfer()` and `cyfer_finalize()`.
#' @param verbose numeric.
#'
#' @returns a list with `bool_converged`, `gene_df` (`gene`, `stat`,
#'   `pvalue`), `lambda` (the selected value) and `z_hat_vec` (the
#'   `cell_imputed_score`, log10 scale, named by cell).
method_cyfer <- function(pca_t1_mat,
                         clone_vec,
                         t2_size_vec,
                         lognorm_t1_mat,
                         lambda_initial = 1,
                         lambda_sequence_length = 10,
                         num_folds = 3,
                         seed_number = 10,
                         verbose = 0){
  stopifnot(nrow(pca_t1_mat) == length(clone_vec),
            !is.null(rownames(pca_t1_mat)),
            setequal(unique(clone_vec), names(t2_size_vec)))
  clone_tab <- table(clone_vec)
  keep_clone_vec <- names(clone_tab)[clone_tab >= 2]
  keep_idx_vec <- which(clone_vec %in% keep_clone_vec)
  x_mat <- scale(pca_t1_mat[keep_idx_vec, , drop = FALSE])
  clone_sub_vec <- as.character(clone_vec[keep_idx_vec])
  lfc_vec <- t2_size_vec[keep_clone_vec]
  num_folds_use <- max(2, min(num_folds, length(keep_clone_vec) - 1))

  fit <- tryCatch({
    fit_res <- multiomeFate::cyfer(cell_features = x_mat,
                                   cell_lineage = clone_sub_vec,
                                   lineage_future_count = lfc_vec,
                                   lambda_initial = lambda_initial,
                                   lambda_sequence_length = lambda_sequence_length,
                                   num_folds = num_folds_use,
                                   seed_number = seed_number,
                                   verbose = max(0, verbose - 1))
    multiomeFate::cyfer_finalize(cell_features = x_mat,
                                 cell_lineage = clone_sub_vec,
                                 fit_res = fit_res,
                                 lineage_future_count = lfc_vec,
                                 seed_number = seed_number)
  }, error = function(e){
    if(verbose > 0) print(paste0("CYFER failed: ", conditionMessage(e)))
    NULL
  })

  num_genes <- ncol(lognorm_t1_mat)
  if(is.null(fit)){
    return(list(bool_converged = FALSE,
                gene_df = data.frame(gene = colnames(lognorm_t1_mat),
                                     stat = rep(NA_real_, num_genes),
                                     pvalue = rep(NA_real_, num_genes),
                                     stringsAsFactors = FALSE),
                lambda = NA_real_,
                z_hat_vec = stats::setNames(rep(NA_real_, nrow(pca_t1_mat)),
                                            rownames(pca_t1_mat))))
  }

  z_hat_vec <- stats::setNames(rep(NA_real_, nrow(pca_t1_mat)),
                               rownames(pca_t1_mat))
  z_hat_vec[keep_idx_vec] <- as.numeric(fit$cell_imputed_score)
  gene_df <- gene_spearman(lognorm_t1_mat[keep_idx_vec, , drop = FALSE],
                           z_hat_vec[keep_idx_vec])

  list(bool_converged = TRUE,
       gene_df = gene_df,
       lambda = fit$lambda,
       z_hat_vec = z_hat_vec)
}

#' Lineage-DE: split clones at the mean t2 size, compare t1 cells of high
#' against low clones gene by gene
#'
#' @param lognorm_t1_mat log-normalized expression, t1 cells by genes.
#' @param clone_vec character vector, the clone of every t1 cell.
#' @param t2_size_vec named numeric vector of t2 clone sizes, zeros kept.
#' @param num_check genes on which the vectorized U is cross-checked
#'   against `stats::wilcox.test()`; 0 skips the check.
#'
#' @returns a list with `gene_df` (`gene`, `stat`, `u`, `pvalue`),
#'   `high_clone_vec`, `num_high_cells` and `num_low_cells`.
method_lineage_de <- function(lognorm_t1_mat,
                              clone_vec,
                              t2_size_vec,
                              num_check = 5){
  stopifnot(nrow(lognorm_t1_mat) == length(clone_vec))
  high_clone_vec <- high_clones(t2_size_vec)
  bool_high_vec <- clone_vec %in% high_clone_vec
  num_genes <- ncol(lognorm_t1_mat)
  if(sum(bool_high_vec) == 0 || sum(!bool_high_vec) == 0){
    return(list(gene_df = data.frame(gene = colnames(lognorm_t1_mat),
                                     stat = rep(NA_real_, num_genes),
                                     u = rep(NA_real_, num_genes),
                                     pvalue = rep(NA_real_, num_genes),
                                     stringsAsFactors = FALSE),
                high_clone_vec = high_clone_vec,
                num_high_cells = sum(bool_high_vec),
                num_low_cells = sum(!bool_high_vec)))
  }
  gene_df <- gene_rank_biserial(lognorm_t1_mat, bool_high_vec)

  if(num_check > 0){
    check_idx_vec <- seq_len(min(num_check, num_genes))
    for(j in check_idx_vec){
      wt <- stats::wilcox.test(lognorm_t1_mat[bool_high_vec, j],
                               lognorm_t1_mat[!bool_high_vec, j],
                               exact = FALSE)
      stopifnot(abs(as.numeric(wt$statistic) - gene_df$u[j]) < 1e-6,
                abs(wt$p.value - gene_df$pvalue[j]) < 1e-6)
    }
  }

  list(gene_df = gene_df,
       high_clone_vec = high_clone_vec,
       num_high_cells = sum(bool_high_vec),
       num_low_cells = sum(!bool_high_vec))
}

#' CoSPAR, matched route: export, run the Python side, correlate every gene
#' with the intraclone fate bias over the t1 cells
#'
#' `state_info` is `"t1"` for t1 cells, `"High"` for t2 cells of the high
#' clones (the lineage-DE split) and `"Low"` for the other t2 cells. The
#' similarity graph is built on the shared PCA. t1 cells of clones with no
#' t2 cells are excluded before the export (memo Section 3.2): they are
#' single-time clones to CoSPAR, outside its transition map, and the 0.5
#' bias it fills in for them would rank the lowest-fate cells in the middle.
#' Their fate bias is returned as `NA` and the gene correlation runs over
#' the included t1 cells only. Fails soft: any error on either side returns
#' `NA` statistics with `bool_success = FALSE`.
#'
#' @param count_mat integer matrix, all cells by genes, with names.
#' @param cell_df data frame with `cell_id`, `time_info`, `clone_id` in the
#'   row order of `count_mat`.
#' @param pca_mat numeric matrix, all cells by PCs, same row order.
#' @param high_clone_vec character vector of the high clones.
#' @param lognorm_t1_mat log-normalized expression, t1 cells by genes.
#' @param export_dir directory for the flat export, under `OUT_ROOT`.
#' @param python_path absolute path to the `cospar` environment's python.
#' @param script_path absolute path to `run_cospar.py`.
#' @param bool_clean remove the bulky export files (counts, PCA, cache)
#'   after a successful import, keeping `manifest.json` and `cospar_out/`.
#' @param data_des cache key for CoSPAR; must differ between datasets.
#' @param seed_number passed to `run_cospar.py --seed`.
#' @param verbose numeric.
#'
#' @returns a list with `bool_success`, `fate_bias_vec` (named over all t1
#'   cells, `NA` for excluded cells and on failure), `gene_df` (`gene`,
#'   `stat`, `pvalue`), `log_vec` (the Python stdout and stderr lines),
#'   `num_high_t2`, `num_low_t2`, `num_progenitor_a`, `num_progenitor_b`,
#'   `num_t1_excluded` and `runtime_sec`.
method_cospar <- function(count_mat,
                          cell_df,
                          pca_mat,
                          high_clone_vec,
                          lognorm_t1_mat,
                          export_dir,
                          python_path,
                          script_path,
                          bool_clean = TRUE,
                          data_des = "dataset",
                          seed_number = 10,
                          verbose = 0){
  stopifnot(nrow(count_mat) == nrow(cell_df), nrow(pca_mat) == nrow(cell_df),
            all(rownames(count_mat) == cell_df$cell_id),
            file.exists(python_path), file.exists(script_path))
  t1_all_id_vec <- cell_df$cell_id[cell_df$time_info == "t1"]

  # Exclude t1 cells of extinct clones (no t2 cells) before the export; see
  # the block comment above. `lognorm_t1_mat` keeps all t1 rows, so the gene
  # step below subsets it by the included cell IDs.
  surviving_clone_vec <- unique(cell_df$clone_id[cell_df$time_info == "t2"])
  bool_keep_vec <- cell_df$time_info == "t2" |
    cell_df$clone_id %in% surviving_clone_vec
  num_t1_excluded <- length(t1_all_id_vec) - sum(bool_keep_vec &
                                                   cell_df$time_info == "t1")
  count_mat <- count_mat[bool_keep_vec, , drop = FALSE]
  pca_mat <- pca_mat[bool_keep_vec, , drop = FALSE]
  cell_df <- cell_df[bool_keep_vec, , drop = FALSE]

  t1_idx <- which(cell_df$time_info == "t1")
  t1_id_vec <- cell_df$cell_id[t1_idx]
  num_genes <- ncol(lognorm_t1_mat)
  empty_gene_df <- data.frame(gene = colnames(lognorm_t1_mat),
                              stat = rep(NA_real_, num_genes),
                              pvalue = rep(NA_real_, num_genes),
                              stringsAsFactors = FALSE)
  na_bias_vec <- stats::setNames(rep(NA_real_, length(t1_all_id_vec)),
                                 t1_all_id_vec)

  state_vec <- ifelse(cell_df$time_info == "t1", "t1",
                      ifelse(cell_df$clone_id %in% high_clone_vec,
                             "High", "Low"))
  num_high_t2 <- sum(state_vec == "High")
  num_low_t2 <- sum(state_vec == "Low")
  fail_list <- list(bool_success = FALSE,
                    fate_bias_vec = na_bias_vec,
                    gene_df = empty_gene_df,
                    log_vec = character(0),
                    num_high_t2 = num_high_t2,
                    num_low_t2 = num_low_t2,
                    num_progenitor_a = NA_integer_,
                    num_progenitor_b = NA_integer_,
                    num_t1_excluded = num_t1_excluded,
                    runtime_sec = NA_real_)
  if(num_high_t2 == 0 || num_low_t2 == 0){
    if(verbose > 0) print("CoSPAR skipped: one t2 fate group is empty")
    return(fail_list)
  }

  start_time <- Sys.time()
  result <- tryCatch({
    export_for_cospar(count_mat = count_mat,
                      time_info = cell_df$time_info,
                      clone_id = cell_df$clone_id,
                      state_info = state_vec,
                      export_dir = export_dir,
                      pca_mat = pca_mat,
                      verbose = max(0, verbose - 1))
    log_vec <- system2(python_path,
                       args = c(shQuote(script_path), shQuote(export_dir),
                                "--t1", "t1", "--t2", "t2",
                                "--fate-a", "High", "--fate-b", "Low",
                                "--data-des", shQuote(data_des),
                                "--seed", seed_number),
                       stdout = TRUE, stderr = TRUE)
    status <- attr(log_vec, "status")
    if(!is.null(status) && status != 0){
      stop("run_cospar.py exited with status ", status, ": ",
           paste(utils::tail(log_vec, 5), collapse = " | "))
    }
    import_list <- import_cospar_results(export_dir)
    list(import_list = import_list, log_vec = log_vec)
  }, error = function(e){
    if(verbose > 0) print(paste0("CoSPAR failed: ", conditionMessage(e)))
    NULL
  })
  runtime_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if(is.null(result)){
    fail_list$runtime_sec <- runtime_sec
    return(fail_list)
  }

  bias_included_vec <- result$import_list$fate_bias[t1_id_vec]
  names(bias_included_vec) <- t1_id_vec
  # A constant bias (every t1 cell at 0.5) is a collapsed map and is reported
  # as missing rather than as a score of 0 (memo Section 7).
  if(anyNA(bias_included_vec) || stats::sd(bias_included_vec) == 0){
    if(verbose > 0) print("CoSPAR failed: fate bias missing or constant")
    fail_list$log_vec <- result$log_vec
    fail_list$runtime_sec <- runtime_sec
    return(fail_list)
  }
  fate_bias_vec <- na_bias_vec
  fate_bias_vec[t1_id_vec] <- bias_included_vec
  gene_df <- gene_spearman(lognorm_t1_mat[t1_id_vec, , drop = FALSE],
                           bias_included_vec)

  if(bool_clean){
    unlink(file.path(export_dir, c("counts.mtx", "X_pca.csv", "cells.txt",
                                   "genes.txt", "obs.csv")))
    unlink(file.path(export_dir, "cospar_cache"), recursive = TRUE)
  }

  list(bool_success = TRUE,
       fate_bias_vec = fate_bias_vec,
       gene_df = gene_df,
       log_vec = result$log_vec,
       num_high_t2 = num_high_t2,
       num_low_t2 = num_low_t2,
       num_progenitor_a = result$import_list$run$n_progenitor_a,
       num_progenitor_b = result$import_list$run$n_progenitor_b,
       num_t1_excluded = num_t1_excluded,
       runtime_sec = runtime_sec)
}
