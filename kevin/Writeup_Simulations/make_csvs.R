# make_csvs.R
# Reads simulation .rds results and writes flat CSVs for downstream use.
# Run: Rscript make_csvs.R

in_dir  <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup_Simulations/"
out_dir <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/csv/kevin/Writeup_Simulations/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(df, name) {
  path <- file.path(out_dir, name)
  write.csv(df, path, row.names = FALSE)
  message("  Wrote: ", name, " (", nrow(df), " rows x ", ncol(df), " cols)")
}

# ==============================================================================
# SIM 1: Growth Modes
# ==============================================================================
message("=== Sim 1: Growth Modes ===")
s1 <- readRDS(file.path(in_dir, "sim1_growth_modes_results.rds"))

# Summary table: one row per growth model (mean/sd/median correlation)
write_csv(s1$summary, "sim1_summary.csv")

# Detailed single-run metrics (correlation, Jaccard, Gini per model)
write_csv(s1$detailed, "sim1_detailed_single_run.csv")

# Replicate-level correlations: one row per replicate, columns = growth models
cor_df <- as.data.frame(s1$replicate_cors)
cor_df$replicate <- seq_len(nrow(cor_df))
cor_df <- cor_df[, c("replicate", setdiff(names(cor_df), "replicate"))]
write_csv(cor_df, "sim1_replicate_cors.csv")

# ==============================================================================
# SIM 2: Rare Resistance
# ==============================================================================
message("=== Sim 2: Rare Resistance ===")
s2 <- readRDS(file.path(in_dir, "sim2_rare_resistance_results.rds"))

# Summary: one row per resistance fraction
write_csv(s2$summary, "sim2_summary.csv")

# Replicate-level details: one row per (resistance_fraction, replicate)
rep_df <- do.call(rbind, lapply(s2$all_results, function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}))
write_csv(rep_df, "sim2_replicate_details.csv")

# ==============================================================================
# SIM 3: Heritability Grid
# ==============================================================================
message("=== Sim 3: Heritability ===")
s3 <- readRDS(file.path(in_dir, "sim3_heritability_results.rds"))

# Summary: one row per (h2_feature, h2_fate) grid cell
write_csv(s3$summary, "sim3_summary.csv")

# ==============================================================================
# SIM 4: Power Analysis
# ==============================================================================
message("=== Sim 4: Power Analysis ===")
s4 <- readRDS(file.path(in_dir, "sim4_power_analysis_results.rds"))

write_csv(s4$axis1_clones,  "sim4_axis1_n_clones.csv")
write_csv(s4$axis2_cellsize, "sim4_axis2_cells_per_clone.csv")
write_csv(s4$axis3_capture,  "sim4_axis3_capture_rate.csv")

# ==============================================================================
# SIM 5: Barcode Dropout
# ==============================================================================
message("=== Sim 5: Barcode Dropout ===")
s5 <- readRDS(file.path(in_dir, "sim5_barcode_dropout_results.rds"))

write_csv(s5$partA, "sim5_partA_clone_dropout.csv")
write_csv(s5$partB, "sim5_partB_cell_subsampling.csv")
if (!is.null(s5$partC)) write_csv(s5$partC, "sim5_partC_gini_vs_clone_size.csv")

# ==============================================================================
# SIM 6: Adaptation Index
# ==============================================================================
message("=== Sim 6: Adaptation Index ===")
s6 <- readRDS(file.path(in_dir, "sim6_adaptation_index_results.rds"))

# Summary: one row per f_dying
write_csv(s6$summary, "sim6_summary.csv")

# Replicate-level details: flatten list-of-lists
rep6 <- do.call(rbind, lapply(s6$all_results, function(rep_list) {
  if (length(rep_list) == 0) return(NULL)
  do.call(rbind, lapply(rep_list, function(r) {
    as.data.frame(r[sapply(r, length) == 1], stringsAsFactors = FALSE)
  }))
}))
write_csv(rep6, "sim6_replicate_details.csv")

# ==============================================================================
# SIM 7: Sensitivity / Specificity
# ==============================================================================
message("=== Sim 7: Sensitivity / Specificity ===")
s7 <- readRDS(file.path(in_dir, "sim7_sensitivity_specificity_results.rds"))

# Summary: one row per scenario (priming vs. plasticity)
write_csv(s7$summary, "sim7_summary.csv")

# Permutation test: CYFER AUROC vs naive AUROC
write_csv(s7$permutation, "sim7_permutation_test.csv")

# Null calibration FPR: one row per null replicate
null_df <- data.frame(
  replicate = seq_along(s7$null_fpr),
  false_positive_rate = s7$null_fpr
)
write_csv(null_df, "sim7_null_calibration_fpr.csv")

# Replicate-level details: flatten nested per-scenario, per-replicate lists
rep7 <- do.call(rbind, lapply(names(s7$all_results), function(sc_name) {
  reps <- s7$all_results[[sc_name]]
  if (length(reps) == 0) return(NULL)
  do.call(rbind, lapply(seq_along(reps), function(i) {
    r <- reps[[i]]
    mc <- r$metrics_cyfer
    mn <- r$metrics_naive
    mcm <- r$metrics_clone_mean

    # Extract clone_mean metrics if available (may be NA scalar)
    is_valid_clone_mean <- is.numeric(mcm) && length(mcm) > 1

    data.frame(
      scenario          = r$scenario,
      replicate         = i,
      seed              = r$seed,
      cor_z_cyfer       = r$cor_z_cyfer,
      auroc_cyfer       = r$auroc_cyfer,
      auroc_naive       = r$auroc_naive,
      auprc_cyfer       = r$auprc_cyfer,
      auprc_naive       = r$auprc_naive,
      # CYFER classification metrics at alpha=0.05
      TP_cyfer          = mc["TP"],
      FP_cyfer          = mc["FP"],
      TN_cyfer          = mc["TN"],
      FN_cyfer          = mc["FN"],
      sensitivity_cyfer = mc["Sensitivity"],
      specificity_cyfer = mc["Specificity"],
      precision_cyfer   = mc["Precision"],
      f1_cyfer          = mc["F1"],
      jaccard_cyfer     = mc["Jaccard"],
      # Naive classification metrics at alpha=0.05
      TP_naive          = mn["TP"],
      FP_naive          = mn["FP"],
      TN_naive          = mn["TN"],
      FN_naive          = mn["FN"],
      sensitivity_naive = mn["Sensitivity"],
      specificity_naive = mn["Specificity"],
      precision_naive   = mn["Precision"],
      f1_naive          = mn["F1"],
      jaccard_naive     = mn["Jaccard"],
      # Clone-mean metrics (NA if unavailable)
      sensitivity_clone_mean = if (is_valid_clone_mean) mcm["Sensitivity"] else NA,
      specificity_clone_mean = if (is_valid_clone_mean) mcm["Specificity"] else NA,
      f1_clone_mean          = if (is_valid_clone_mean) mcm["F1"]          else NA,
      jaccard_clone_mean     = if (is_valid_clone_mean) mcm["Jaccard"]     else NA,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}))
write_csv(rep7, "sim7_replicate_details.csv")

message("\nDone. All CSVs written to:\n  ", out_dir)
