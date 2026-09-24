# func_sweep_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-22
#
# The sweep engine shared by `sim_gini_claude.R` and `sim_aed_claude.R`
# (simulation_design_claude.md, Sections 5, 7 and 8): calibrate every level,
# then for each level x replicate generate a dataset, run the three methods,
# score them, and append a time-stamped line to a progress file at every
# stage. The RDS is re-saved after every level so a killed run keeps what it
# finished.
#
# Requires `func_generate_claude.R`, `func_methods_claude.R` and
# `cospar_flat_io.R` to be sourced first.

# Progress reporting ----------------------------------------------------------

#' Append a time-stamped line to the progress file and print it
#'
#' @param message_str the line.
#' @param progress_file path; `NULL` prints only.
#' @param start_time the run's start, for the elapsed prefix.
.log_progress <- function(message_str, progress_file, start_time){
  elapsed_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  line_str <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                     sprintf("%6.1f min", elapsed_min), " | ", message_str)
  print(line_str)
  if(!is.null(progress_file)){
    cat(line_str, "\n", file = progress_file, append = TRUE, sep = "")
  }
  invisible(line_str)
}

# One dataset -----------------------------------------------------------------

#' Generate one dataset, run the three methods, score them
#'
#' @param h2,latent_scale,tau_delta,spread_variation,kappa generator knobs.
#' @param seed_number seed for the generator, CYFER and CoSPAR.
#' @param export_dir directory for this dataset's CoSPAR export.
#' @param python_path,script_path for `method_cospar()`.
#' @param d_latent latent dimension of the generator.
#' @param d_pca number of PCs in the shared embedding.
#' @param data_des CoSPAR cache key, unique per dataset.
#' @param generator_args list of further arguments to `generate_dataset()`.
#' @param progress_file,start_time for `.log_progress()`.
#' @param stage_prefix string prefixed to every progress line.
#' @param verbose numeric.
#'
#' @returns a list with `detail_df` (three rows, one per method, carrying the
#'   metrics, score-level diagnostics and every dataset-level realized
#'   statistic), `gene_df` (per-gene statistics of the truth and the three
#'   methods) and `cospar_log_vec`.
run_one_dataset <- function(h2,
                            latent_scale,
                            tau_delta,
                            spread_variation,
                            kappa,
                            seed_number,
                            export_dir,
                            python_path,
                            script_path,
                            d_latent = 10,
                            d_pca = 10,
                            data_des = "dataset",
                            generator_args = list(),
                            progress_file = NULL,
                            start_time = Sys.time(),
                            stage_prefix = "",
                            verbose = 0){
  timer_vec <- c(generate = NA_real_, embed = NA_real_, cyfer = NA_real_,
                 lineage_de = NA_real_, cospar = NA_real_)

  # Generation and embedding.
  t0 <- Sys.time()
  dat <- do.call(generate_dataset,
                 c(list(h2 = h2, latent_scale = latent_scale, d = d_latent,
                        kappa = kappa, spread_variation = spread_variation,
                        tau_delta = tau_delta, seed_number = seed_number),
                   generator_args))
  timer_vec["generate"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  emb <- compute_embedding(dat$count_mat, dat$t1_idx, d = d_pca)
  pca_t1_mat <- emb$pca_mat[dat$t1_idx, , drop = FALSE]
  lognorm_t1_mat <- emb$lognorm_mat[dat$t1_idx, , drop = FALSE]
  clone_vec <- dat$cell_df$clone_id[dat$t1_idx]
  t2_size_vec <- dat$t2_size_vec
  aed_vec <- compute_aed(pca_t1_mat, clone_vec)
  truth_list <- compute_truth(lognorm_t1_mat, dat$z_true_vec)
  timer_vec["embed"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  .log_progress(paste0(stage_prefix, "generated + embedded: t2 cells = ",
                       sum(t2_size_vec), ", Gini = ",
                       round(gini_coef(t2_size_vec), 3), ", mean AED = ",
                       round(mean(aed_vec), 3), " (",
                       round(timer_vec["generate"] + timer_vec["embed"], 1),
                       " s)"),
                progress_file, start_time)

  # CYFER.
  t0 <- Sys.time()
  cyfer_list <- method_cyfer(pca_t1_mat, clone_vec, t2_size_vec,
                             lognorm_t1_mat, seed_number = seed_number,
                             verbose = verbose)
  timer_vec["cyfer"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  .log_progress(paste0(stage_prefix, "CYFER done: converged = ",
                       cyfer_list$bool_converged, ", lambda = ",
                       signif(cyfer_list$lambda, 3), " (",
                       round(timer_vec["cyfer"], 1), " s)"),
                progress_file, start_time)

  # Lineage-DE.
  t0 <- Sys.time()
  de_list <- method_lineage_de(lognorm_t1_mat, clone_vec, t2_size_vec)
  timer_vec["lineage_de"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # CoSPAR.
  t0 <- Sys.time()
  cospar_list <- method_cospar(count_mat = dat$count_mat,
                               cell_df = dat$cell_df,
                               pca_mat = emb$pca_mat,
                               high_clone_vec = de_list$high_clone_vec,
                               lognorm_t1_mat = lognorm_t1_mat,
                               export_dir = export_dir,
                               python_path = python_path,
                               script_path = script_path,
                               data_des = data_des,
                               seed_number = seed_number,
                               verbose = verbose)
  timer_vec["cospar"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  .log_progress(paste0(stage_prefix, "CoSPAR done: success = ",
                       cospar_list$bool_success, ", High/Low t2 = ",
                       cospar_list$num_high_t2, "/", cospar_list$num_low_t2,
                       " (", round(timer_vec["cospar"], 1), " s)"),
                progress_file, start_time)

  # Scores.
  z_true_vec <- dat$z_true_vec
  score_row <- function(method, gene_df, score_vec, bool_converged,
                        runtime_sec){
    metric_df <- score_gene_stats(gene_df$stat, gene_df$pvalue, truth_list)
    bool_score <- !all(is.na(score_vec))
    data.frame(method = method,
               metric_df,
               score_cor = if(bool_score) stats::cor(score_vec, z_true_vec,
                                                     use = "complete.obs") else NA_real_,
               score_spearman = if(bool_score) stats::cor(score_vec, z_true_vec,
                                                          method = "spearman",
                                                          use = "complete.obs") else NA_real_,
               bool_converged = bool_converged,
               runtime_sec = runtime_sec,
               stringsAsFactors = FALSE)
  }
  detail_df <- rbind(
    score_row("CYFER", cyfer_list$gene_df, cyfer_list$z_hat_vec,
              cyfer_list$bool_converged, timer_vec["cyfer"]),
    score_row("CoSPAR", cospar_list$gene_df, cospar_list$fate_bias_vec,
              cospar_list$bool_success, timer_vec["cospar"]),
    score_row("Lineage-DE", de_list$gene_df, NA_real_,
              !all(is.na(de_list$gene_df$stat)), timer_vec["lineage_de"]))

  # Dataset-level realized statistics and diagnostics (Sections 5 and 7).
  dataset_df <- data.frame(
    seed = seed_number,
    h2 = h2,
    latent_scale = latent_scale,
    beta_0 = dat$latent$beta_0,
    gini_t2 = gini_coef(t2_size_vec),
    gini_pooled = gini_coef(t2_size_vec + dat$param_list$num_cells_per_clone),
    aed_mean = mean(aed_vec),
    aed_min = min(aed_vec),
    aed_q05 = unname(stats::quantile(aed_vec, 0.05)),
    aed_median = stats::median(aed_vec),
    aed_q95 = unname(stats::quantile(aed_vec, 0.95)),
    aed_max = max(aed_vec),
    t2_total = sum(t2_size_vec),
    t2_max_clone = max(t2_size_vec),
    t2_num_zero = sum(t2_size_vec == 0),
    spearman_aed_t2size = stats::cor(aed_vec[names(t2_size_vec)], t2_size_vec,
                                     method = "spearman"),
    pca_heritability = compute_pca_heritability(pca_t1_mat, clone_vec),
    num_high_clones = length(de_list$high_clone_vec),
    num_high_cells = de_list$num_high_cells,
    num_low_cells = de_list$num_low_cells,
    num_high_t2 = cospar_list$num_high_t2,
    num_low_t2 = cospar_list$num_low_t2,
    num_progenitor_a = cospar_list$num_progenitor_a,
    num_progenitor_b = cospar_list$num_progenitor_b,
    num_t1_excluded_cospar = cospar_list$num_t1_excluded,
    time_generate_sec = timer_vec["generate"],
    time_embed_sec = timer_vec["embed"],
    stringsAsFactors = FALSE)
  detail_df <- cbind(detail_df, dataset_df[rep(1, nrow(detail_df)), ],
                     row.names = NULL)

  gene_df <- data.frame(gene = truth_list$truth_df$gene,
                        truth = truth_list$truth_df$stat,
                        truth_qvalue = truth_list$truth_df$qvalue,
                        cyfer = cyfer_list$gene_df$stat,
                        cyfer_pvalue = cyfer_list$gene_df$pvalue,
                        cospar = cospar_list$gene_df$stat,
                        cospar_pvalue = cospar_list$gene_df$pvalue,
                        lineage_de = de_list$gene_df$stat,
                        lineage_de_pvalue = de_list$gene_df$pvalue,
                        stringsAsFactors = FALSE)

  list(cospar_log_vec = cospar_list$log_vec,
       detail_df = detail_df,
       gene_df = gene_df)
}

# The sweep -------------------------------------------------------------------

#' Run one axis of the Figure 4 sweep
#'
#' @param axis `"gini"` (levels are t2 Gini targets, AED held) or `"aed"`
#'   (levels are mean squared AED targets, Gini held).
#' @param level_vec the targets, in order.
#' @param fixed_val the held statistic: the mean AED for `"gini"`, the t2
#'   Gini for `"aed"`.
#' @param num_replicates replicates per level.
#' @param out_dir directory for the RDS, the CoSPAR exports and the progress
#'   file.
#' @param python_path,script_path for `method_cospar()`.
#' @param rds_file path of the RDS to write.
#' @param progress_file path of the progress text file.
#' @param aed_tol,max_recentre the realized mean AED of a check dataset is
#'   compared with the level's AED target and `h2` shifted by the miss up to
#'   `max_recentre` times when it exceeds `aed_tol` (memo Section 5).
#' @param d_pca number of PCs in the shared embedding (the generator's
#'   latent dimension stays at its default of 10).
#' @param generator_args further arguments to `generate_dataset()`.
#' @param kappa,spread_variation,tau_delta generator knobs fixed across the
#'   sweep.
#' @param run_label names this run's CoSPAR export directories and cache
#'   keys, so two runs of the same axis with different settings do not
#'   overwrite each other's exports; defaults to `axis`.
#' @param seed_base replicate `r` of level `j` uses seed
#'   `seed_base + 100 * j + r`.
#' @param verbose numeric.
#'
#' @returns the list saved to `rds_file`: `summary` (one row per level x
#'   method), `replicate_details` (one row per level x replicate x method),
#'   `calibration` (one row per level), `gene_stats` (list of per-dataset
#'   gene data frames, named `level<j>_rep<r>`), `cospar_logs` and `params`.
run_sweep <- function(axis,
                      level_vec,
                      fixed_val,
                      num_replicates,
                      out_dir,
                      python_path,
                      script_path,
                      rds_file,
                      progress_file,
                      aed_tol = 0.03,
                      d_pca = 10,
                      generator_args = list(),
                      max_recentre = 3,
                      kappa = 1.5,
                      run_label = axis,
                      spread_variation = "noncausal",
                      tau_delta = 0.6,
                      seed_base = 0,
                      verbose = 0){
  stopifnot(axis %in% c("gini", "aed"), length(level_vec) > 0,
            num_replicates >= 1)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  start_time <- Sys.time()
  cat("", file = progress_file)
  .log_progress(paste0("=== ", axis, " sweep: ", length(level_vec),
                       " levels x ", num_replicates, " replicates = ",
                       length(level_vec) * num_replicates, " datasets, fixed ",
                       if(axis == "gini") "mean AED = " else "t2 Gini = ",
                       fixed_val, ", tau_delta = ", tau_delta,
                       ", spread_variation = ", spread_variation,
                       ", kappa = ", kappa, ", d_pca = ", d_pca),
                progress_file, start_time)
  .log_progress(paste0("progress file: ", progress_file), progress_file,
                start_time)
  .log_progress(paste0("RDS (re-saved after every level): ", rds_file),
                progress_file, start_time)

  # Calibration (Section 5): h2 from the AED target, s bisected for the
  # Gini target, both per level.
  calibration_list <- vector("list", length(level_vec))
  for(j in seq_along(level_vec)){
    if(axis == "gini"){
      target_gini <- level_vec[j]
      target_aed <- fixed_val
    } else {
      target_gini <- fixed_val
      target_aed <- level_vec[j]
    }
    # Calibration draws use seeds seed_base + 100 j + 50 + k, disjoint from
    # the replicates' seed_base + 100 j + r (r at most 20).
    calib_seed <- seed_base + 100 * j + 50
    h2_initial <- 1 - target_aed
    h2 <- h2_initial
    aed_check <- NA_real_
    num_recentre <- 0
    # Section 5: the mean AED is about 1 - h2 plus a count-noise floor that
    # grows as the latent scale shrinks (and with the number of PCs); check
    # the realized mean on one full dataset and shift h2 by the miss, up to
    # max_recentre times, re-bisecting s each time. h2 is clamped to [0, 1],
    # so a target below the floor is reported as missed, not silently moved.
    for(recentre in seq_len(max_recentre + 1)){
      calib <- calibrate_scale(h2 = h2,
                               gini_target = target_gini,
                               kappa = kappa,
                               spread_variation = spread_variation,
                               seed_number = calib_seed)
      check_dat <- do.call(generate_dataset,
                           c(list(h2 = h2, latent_scale = calib$latent_scale,
                                  kappa = kappa,
                                  spread_variation = spread_variation,
                                  tau_delta = tau_delta,
                                  seed_number = calib_seed),
                             generator_args))
      check_emb <- compute_embedding(check_dat$count_mat, check_dat$t1_idx,
                                     d = d_pca)
      aed_check <- mean(compute_aed(
        check_emb$pca_mat[check_dat$t1_idx, , drop = FALSE],
        check_dat$cell_df$clone_id[check_dat$t1_idx]))
      if(abs(aed_check - target_aed) <= aed_tol || recentre > max_recentre) break
      num_recentre <- recentre
      h2 <- min(1, max(0, h2 + (aed_check - target_aed)))
    }
    calibration_list[[j]] <- data.frame(level_idx = j,
                                        level = level_vec[j],
                                        target_gini = target_gini,
                                        target_aed = target_aed,
                                        h2_initial = h2_initial,
                                        h2 = h2,
                                        num_recentre = num_recentre,
                                        aed_check = aed_check,
                                        latent_scale = calib$latent_scale,
                                        tau = calib$tau,
                                        sigma_w = calib$sigma_w,
                                        gini_calibrated = calib$gini_realized,
                                        num_iter = calib$num_iter)
    .log_progress(paste0("calibrated level ", j, " (", level_vec[j],
                         "): h2 = ", round(h2, 3), " (", num_recentre,
                         " re-centring steps from ", round(h2_initial, 3),
                         "), s = ", signif(calib$latent_scale, 4), ", tau = ",
                         signif(calib$tau, 3), ", sigma_w = ",
                         signif(calib$sigma_w, 3), ", mean Gini over draws = ",
                         round(calib$gini_realized, 3), ", check AED = ",
                         round(aed_check, 3)),
                  progress_file, start_time)
  }
  calibration_df <- do.call(rbind, calibration_list)

  params <- list(aed_tol = aed_tol, axis = axis, d_pca = d_pca,
                 fixed_val = fixed_val, max_recentre = max_recentre,
                 generator_args = generator_args, kappa = kappa,
                 level_vec = level_vec, num_replicates = num_replicates,
                 python_path = python_path, run_label = run_label,
                 script_path = script_path, seed_base = seed_base, spread_variation = spread_variation,
                 tau_delta = tau_delta,
                 multiomeFate_version = as.character(
                   utils::packageVersion("multiomeFate")),
                 r_version = R.version.string,
                 start_time = start_time)

  # Level x replicate loop, saving after every level.
  num_datasets <- length(level_vec) * num_replicates
  detail_list <- list()
  gene_list <- list()
  cospar_log_list <- list()
  duration_vec <- numeric(0)
  for(j in seq_along(level_vec)){
    for(r in seq_len(num_replicates)){
      dataset_name <- paste0("level", j, "_rep", r)
      seed_number <- seed_base + 100 * j + r
      prefix <- paste0("[", axis, " level ", j, "/", length(level_vec),
                       " (", level_vec[j], "), rep ", r, "/", num_replicates,
                       "] ")
      .log_progress(paste0(prefix, "start, seed = ", seed_number),
                    progress_file, start_time)
      t0 <- Sys.time()

      one_list <- run_one_dataset(
        h2 = calibration_df$h2[j],
        latent_scale = calibration_df$latent_scale[j],
        tau_delta = tau_delta,
        spread_variation = spread_variation,
        kappa = kappa,
        seed_number = seed_number,
        d_pca = d_pca,
        export_dir = file.path(out_dir, "cospar_exports",
                               paste0(run_label, "_", dataset_name)),
        python_path = python_path,
        script_path = script_path,
        data_des = paste0(run_label, "_", dataset_name),
        generator_args = generator_args,
        progress_file = progress_file,
        start_time = start_time,
        stage_prefix = prefix,
        verbose = verbose)

      detail_df <- one_list$detail_df
      detail_df <- cbind(data.frame(level_idx = j, level = level_vec[j],
                                    replicate = r,
                                    target_gini = calibration_df$target_gini[j],
                                    target_aed = calibration_df$target_aed[j]),
                         detail_df, row.names = NULL)
      detail_list[[dataset_name]] <- detail_df
      gene_list[[dataset_name]] <- one_list$gene_df
      cospar_log_list[[dataset_name]] <- one_list$cospar_log_vec

      duration_vec <- c(duration_vec,
                        as.numeric(difftime(Sys.time(), t0, units = "mins")))
      num_done <- length(duration_vec)
      remaining_min <- mean(duration_vec) * (num_datasets - num_done)
      metric_str <- paste0(detail_df$method, " = ",
                           round(detail_df$spearman, 3), collapse = ", ")
      .log_progress(paste0(prefix, "done in ",
                           round(duration_vec[num_done], 2), " min; Spearman: ",
                           metric_str, "; ", num_done, "/", num_datasets,
                           " datasets, about ", round(remaining_min, 1),
                           " min remaining"),
                    progress_file, start_time)
    }

    replicate_details_df <- do.call(rbind, detail_list)
    rownames(replicate_details_df) <- NULL
    summary_df <- .summarize_sweep(replicate_details_df)
    result_list <- list(calibration = calibration_df,
                        cospar_logs = cospar_log_list,
                        gene_stats = gene_list,
                        params = params,
                        replicate_details = replicate_details_df,
                        summary = summary_df)
    saveRDS(result_list, rds_file)
    # A compact per-level line, so the curve can be read off the progress
    # file with `tail -f` during a long run without loading the RDS.
    level_summary_df <- summary_df[summary_df$level_idx == j, ]
    level_str <- paste0(level_summary_df$method, " = ",
                        round(level_summary_df$spearman_mean, 3), " (sd ",
                        round(level_summary_df$spearman_sd, 3), ")",
                        collapse = ", ")
    .log_progress(paste0("*** level ", j, "/", length(level_vec), " (",
                         level_vec[j], ") complete over ",
                         level_summary_df$num_replicates[1],
                         " replicates; mean Spearman: ", level_str),
                  progress_file, start_time)
    .log_progress(paste0("saved ", rds_file, " through level ", j),
                  progress_file, start_time)
  }

  .log_progress(paste0("=== ", axis, " sweep finished"), progress_file,
                start_time)
  result_list
}

#' Mean and SD of every metric per level x method
#'
#' @param replicate_details_df the per-replicate table of `run_sweep()`.
#'
#' @returns a data frame, one row per level x method.
#' @noRd
.summarize_sweep <- function(replicate_details_df){
  key_df <- unique(replicate_details_df[, c("level_idx", "level", "method",
                                            "target_gini", "target_aed")])
  key_df <- key_df[order(key_df$level_idx, key_df$method), ]
  metric_vec <- c("spearman", "pearson", "jaccard", "score_cor",
                  "score_spearman")
  realized_vec <- c("gini_t2", "gini_pooled", "aed_mean", "aed_q05",
                    "aed_q95", "t2_total", "t2_max_clone", "t2_num_zero",
                    "spearman_aed_t2size", "pca_heritability",
                    "num_high_clones", "num_t1_excluded_cospar",
                    "runtime_sec")
  summary_list <- lapply(seq_len(nrow(key_df)), function(i){
    sub_df <- replicate_details_df[
      replicate_details_df$level_idx == key_df$level_idx[i] &
        replicate_details_df$method == key_df$method[i], ]
    row_df <- key_df[i, ]
    row_df$num_replicates <- nrow(sub_df)
    row_df$num_converged <- sum(sub_df$bool_converged)
    for(metric in metric_vec){
      row_df[[paste0(metric, "_mean")]] <- mean(sub_df[[metric]], na.rm = TRUE)
      row_df[[paste0(metric, "_sd")]] <- stats::sd(sub_df[[metric]], na.rm = TRUE)
    }
    for(realized in realized_vec){
      row_df[[paste0(realized, "_mean")]] <- mean(sub_df[[realized]],
                                                  na.rm = TRUE)
    }
    row_df
  })
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  summary_df
}
