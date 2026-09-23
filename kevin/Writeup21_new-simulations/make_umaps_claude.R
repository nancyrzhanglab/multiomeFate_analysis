# make_umaps_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-23
#
# UMAPs of example datasets from the two Writeup21 sweeps, for the
# trailblazing report: log-normalize, PCA (the shared 20-PC embedding of
# simulation_design_claude.md, Section 6.2), then UMAP, for the easiest,
# middle and hardest level of each axis (Gini 0.3/0.6/0.9, AED 0.1/0.55/1.0),
# replicate 1 of each. The datasets are regenerated exactly from the `_r2`
# calibration tables and the sweep's seed scheme, so `sim_gini_claude.R` and
# `sim_aed_claude.R` must have run first.
#
# Run:  Rscript make_umaps_claude.R
# Writes csv/kevin/Writeup21/umap_coords_r2.csv, one row per cell per
# dataset, which the report colors by time point, clone and true fate
# potential (t1 cells only).

library(uwot)

rm(list = ls())

# Paths ------------------------------------------------------------------------

repo_dir <- file.path("/Users/kevinlin/Library/CloudStorage/Dropbox",
                      "Collaboration-and-People/archive/Nancy/multiomeFate",
                      "git/multiomeFate_analysis")
out_dir <- file.path("/Users/kevinlin/Library/CloudStorage/Dropbox",
                     "Collaboration-and-People/archive/Nancy/multiomeFate",
                     "out/Writeup21_new-simulations")
script_dir <- file.path(repo_dir, "kevin", "Writeup21_new-simulations")
csv_dir <- file.path(repo_dir, "csv", "kevin", "Writeup21")

source(file.path(script_dir, "func_generate_claude.R"))
source(file.path(script_dir, "func_methods_claude.R"))

# Settings ---------------------------------------------------------------------

d_pca <- 20
tau_delta <- 0.6
kappa <- 1.5
spread_variation <- "noncausal"
# Easiest / middle / hardest level of each seven-level ladder, replicate 1;
# each axis's seed_base matches its driver.
setting_df <- data.frame(axis = rep(c("gini", "aed"), each = 3),
                         level_idx = rep(c(1, 4, 7), times = 2),
                         difficulty = rep(c("easiest", "middle", "hardest"),
                                          times = 2),
                         seed_base = rep(c(1000, 2000), each = 3),
                         stringsAsFactors = FALSE)

# Run ---------------------------------------------------------------------------

umap_list <- list()
for(i in seq_len(nrow(setting_df))){
  axis <- setting_df$axis[i]
  level_idx <- setting_df$level_idx[i]
  rds_file <- file.path(out_dir, paste0("sim_", axis, "_claude_r2.rds"))
  stopifnot(file.exists(rds_file))
  calibration_df <- readRDS(rds_file)$calibration
  level <- calibration_df$level[level_idx]
  seed_number <- setting_df$seed_base[i] + 100 * level_idx + 1

  print(paste0("Axis ", axis, ", level ", level, " (",
               setting_df$difficulty[i], "), seed ", seed_number))
  dat <- generate_dataset(h2 = calibration_df$h2[level_idx],
                          latent_scale = calibration_df$latent_scale[level_idx],
                          kappa = kappa,
                          spread_variation = spread_variation,
                          tau_delta = tau_delta,
                          seed_number = seed_number)
  emb <- compute_embedding(dat$count_mat, dat$t1_idx, d = d_pca)

  set.seed(seed_number)
  umap_mat <- uwot::umap(emb$pca_mat)

  z_true_vec <- rep(NA_real_, nrow(dat$cell_df))
  z_true_vec[dat$t1_idx] <- dat$z_true_vec
  umap_list[[i]] <- data.frame(axis = axis,
                               level_idx = level_idx,
                               level = level,
                               difficulty = setting_df$difficulty[i],
                               cell_id = dat$cell_df$cell_id,
                               time_info = dat$cell_df$time_info,
                               clone_id = dat$cell_df$clone_id,
                               z_true = z_true_vec,
                               umap1 = umap_mat[, 1],
                               umap2 = umap_mat[, 2],
                               stringsAsFactors = FALSE)
}

# Export -----------------------------------------------------------------------

umap_df <- do.call(rbind, umap_list)
rownames(umap_df) <- NULL
umap_df$z_true <- signif(umap_df$z_true, 5)
umap_df$umap1 <- signif(umap_df$umap1, 5)
umap_df$umap2 <- signif(umap_df$umap2, 5)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(umap_df, file.path(csv_dir, "umap_coords_r2.csv"),
                 row.names = FALSE)
print(paste0("Wrote ", file.path(csv_dir, "umap_coords_r2.csv"), " (",
             nrow(umap_df), " rows)"))
