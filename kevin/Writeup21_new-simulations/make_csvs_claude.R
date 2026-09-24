# make_csvs_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-22
#
# Reads the two sweep RDS files written by `sim_gini_claude.R` and
# `sim_aed_claude.R` from OUT_ROOT and writes flat CSVs to
# csv/kevin/Writeup21/. A missing RDS is skipped with a message.
#
# Run:   Rscript make_csvs_claude.R                     (reads the `_r2` RDS)
# Round: WRITEUP21_ROUND=final Rscript make_csvs_claude.R  (reads `_final`,
#        the 20-replicate round of `sim_final_claude.R`)
# Dry:   WRITEUP21_DRY=1 Rscript make_csvs_claude.R    (reads the `_dry` RDS)
# PCs:   WRITEUP21_PCA=<d> Rscript make_csvs_claude.R  (reads `_pca<d>`)
#
# Per axis:
#   sim_<axis>_summary.csv             one row per level x method, mean/SD
#   sim_<axis>_replicate_details.csv   one row per level x replicate x method
#   sim_<axis>_calibration.csv         one row per level

rm(list = ls())

# Paths ------------------------------------------------------------------------

repo_dir <- file.path("/Users/kevinlin/Library/CloudStorage/Dropbox",
                      "Collaboration-and-People/archive/Nancy/multiomeFate",
                      "git/multiomeFate_analysis")
out_dir <- file.path("/Users/kevinlin/Library/CloudStorage/Dropbox",
                     "Collaboration-and-People/archive/Nancy/multiomeFate",
                     "out/Writeup21_new-simulations")
csv_dir <- file.path(repo_dir, "csv", "kevin", "Writeup21")

bool_dry <- Sys.getenv("WRITEUP21_DRY") == "1"
d_pca <- as.numeric(Sys.getenv("WRITEUP21_PCA", unset = "20"))
# `r2` is the trailblazing round (2 replicates), `final` the 20-replicate
# round; both stay on disk, so the round names the files rather than
# replacing them.
round_str <- Sys.getenv("WRITEUP21_ROUND", unset = "r2")
suffix <- paste0("_", round_str,
                 if(d_pca != 20) paste0("_pca", d_pca) else "",
                 if(bool_dry) "_dry" else "")

# Export -----------------------------------------------------------------------

dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)

for(axis in c("gini", "aed")){
  rds_file <- file.path(out_dir, paste0("sim_", axis, "_claude", suffix,
                                        ".rds"))
  if(!file.exists(rds_file)){
    print(paste0("Skipping ", axis, ": ", rds_file, " not found"))
    next
  }
  result_list <- readRDS(rds_file)

  summary_df <- result_list$summary
  numeric_col_vec <- sapply(summary_df, is.numeric)
  summary_df[numeric_col_vec] <- lapply(summary_df[numeric_col_vec],
                                        signif, digits = 4)
  utils::write.csv(summary_df,
                   file.path(csv_dir, paste0("sim_", axis, "_summary",
                                             suffix, ".csv")),
                   row.names = FALSE)

  detail_df <- result_list$replicate_details
  numeric_col_vec <- sapply(detail_df, is.numeric)
  detail_df[numeric_col_vec] <- lapply(detail_df[numeric_col_vec],
                                       signif, digits = 4)
  utils::write.csv(detail_df,
                   file.path(csv_dir, paste0("sim_", axis,
                                             "_replicate_details", suffix,
                                             ".csv")),
                   row.names = FALSE)

  utils::write.csv(result_list$calibration,
                   file.path(csv_dir, paste0("sim_", axis, "_calibration",
                                             suffix, ".csv")),
                   row.names = FALSE)

  print(paste0("Wrote ", axis, " CSVs (", nrow(summary_df), " summary rows, ",
               nrow(detail_df), " replicate rows)"))
}
