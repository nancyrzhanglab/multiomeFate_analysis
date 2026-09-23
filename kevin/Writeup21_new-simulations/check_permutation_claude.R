# check_permutation_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-23
#
# The label-permutation sanity check of simulation_design_claude.md, Section 7:
# at the middle level of each axis (Gini 0.6, AED 0.55), permute the t1
# cells' clone labels before fitting, so no feature carries fate information,
# and confirm that every method's metric sits near 0. The t2 clone sizes and
# the t2 cells' own clone labels are untouched: the permutation severs the
# t1-side link between expression and fate, which is exactly what each method
# relies on (CYFER's clone counts, lineage-DE's split, CoSPAR's barcode link).
#
# Reads the calibrated (h2, s) for the middle levels from the `_r2` sweep RDS
# files, so `sim_gini_claude.R` and `sim_aed_claude.R` must have run first.
#
# Run:  Rscript check_permutation_claude.R
# Writes csv/kevin/Writeup21/check_permutation_r2.csv (~1 min per dataset,
# 2 axes x 2 seeds).

library(multiomeFate)
library(Matrix)
library(jsonlite)

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
python_path <- "/opt/miniconda3/envs/cospar/bin/python"
script_path <- file.path(script_dir, "run_cospar.py")

source(file.path(script_dir, "func_generate_claude.R"))
source(file.path(script_dir, "func_methods_claude.R"))
source(file.path(script_dir, "func_sweep_claude.R"))
source(file.path(script_dir, "cospar_flat_io.R"))

# Settings ---------------------------------------------------------------------

d_pca <- 20
# Five permutations per axis rather than one or two: under the null CYFER's
# cross-validation does not shrink the fit all the way to a constant, so its
# permuted score is an overfit *random* direction in the 20-PC space, whose
# chance alignment with Z_true has a standard deviation near 1/sqrt(20); a
# handful of permutations is needed to see that spread rather than one draw.
num_seeds <- 5
tau_delta <- 0.6
kappa <- 1.5
spread_variation <- "noncausal"
# The middle level of each seven-level ladder (index 4), and each axis's
# seed_base from its driver. The permutation seeds sit at offset 80 + k,
# disjoint from the replicates (+1..+20) and the calibration draws (+50 + k).
axis_df <- data.frame(axis = c("gini", "aed"),
                      level_idx = c(4, 4),
                      seed_base = c(1000, 2000),
                      stringsAsFactors = FALSE)

# Run ---------------------------------------------------------------------------

result_list <- list()
for(i in seq_len(nrow(axis_df))){
  axis <- axis_df$axis[i]
  level_idx <- axis_df$level_idx[i]
  rds_file <- file.path(out_dir, paste0("sim_", axis, "_claude_r2.rds"))
  stopifnot(file.exists(rds_file))
  calibration_df <- readRDS(rds_file)$calibration
  h2 <- calibration_df$h2[level_idx]
  latent_scale <- calibration_df$latent_scale[level_idx]
  level <- calibration_df$level[level_idx]
  print(paste0("Axis ", axis, ", level ", level, ": h2 = ", round(h2, 3),
               ", s = ", signif(latent_scale, 4)))

  for(k in seq_len(num_seeds)){
    seed_number <- axis_df$seed_base[i] + 100 * level_idx + 80 + k
    dat <- generate_dataset(h2 = h2,
                            latent_scale = latent_scale,
                            kappa = kappa,
                            spread_variation = spread_variation,
                            tau_delta = tau_delta,
                            seed_number = seed_number)
    emb <- compute_embedding(dat$count_mat, dat$t1_idx, d = d_pca)
    pca_t1_mat <- emb$pca_mat[dat$t1_idx, , drop = FALSE]
    lognorm_t1_mat <- emb$lognorm_mat[dat$t1_idx, , drop = FALSE]
    truth_list <- compute_truth(lognorm_t1_mat, dat$z_true_vec)

    # The permutation: reassign the t1 cells to clones uniformly at random.
    # The t2 sizes stay attached to the clone names, so every method sees the
    # same marginal clone-size distribution with a severed feature-fate link.
    set.seed(seed_number)
    clone_perm_vec <- sample(dat$cell_df$clone_id[dat$t1_idx])
    cell_perm_df <- dat$cell_df
    cell_perm_df$clone_id[dat$t1_idx] <- clone_perm_vec

    cyfer_list <- method_cyfer(pca_t1_mat, clone_perm_vec, dat$t2_size_vec,
                               lognorm_t1_mat, seed_number = seed_number)
    de_list <- method_lineage_de(lognorm_t1_mat, clone_perm_vec,
                                 dat$t2_size_vec)
    cospar_list <- method_cospar(count_mat = dat$count_mat,
                                 cell_df = cell_perm_df,
                                 pca_mat = emb$pca_mat,
                                 high_clone_vec = de_list$high_clone_vec,
                                 lognorm_t1_mat = lognorm_t1_mat,
                                 export_dir = file.path(
                                   out_dir, "cospar_exports",
                                   paste0("perm_", axis, "_seed", k)),
                                 python_path = python_path,
                                 script_path = script_path,
                                 data_des = paste0("perm_", axis, "_seed", k),
                                 seed_number = seed_number)

    score_one <- function(method, gene_df, score_vec, bool_converged){
      metric_df <- score_gene_stats(gene_df$stat, gene_df$pvalue, truth_list)
      bool_score <- !all(is.na(score_vec))
      data.frame(axis = axis,
                 level = level,
                 seed = seed_number,
                 method = method,
                 metric_df[, c("spearman", "pearson", "jaccard")],
                 score_cor = if(bool_score){
                   stats::cor(score_vec, dat$z_true_vec, use = "complete.obs")
                 } else NA_real_,
                 bool_converged = bool_converged,
                 stringsAsFactors = FALSE)
    }
    row_df <- rbind(
      score_one("CYFER", cyfer_list$gene_df, cyfer_list$z_hat_vec,
                cyfer_list$bool_converged),
      score_one("CoSPAR", cospar_list$gene_df, cospar_list$fate_bias_vec,
                cospar_list$bool_success),
      score_one("Lineage-DE", de_list$gene_df, NA_real_,
                !all(is.na(de_list$gene_df$stat))))
    result_list[[paste0(axis, "_seed", k)]] <- row_df
    print(paste0("  seed ", seed_number, ": Spearman ",
                 paste0(row_df$method, " = ", round(row_df$spearman, 3),
                        collapse = ", ")))
  }
}

# Export -----------------------------------------------------------------------

check_df <- do.call(rbind, result_list)
rownames(check_df) <- NULL
numeric_col_vec <- sapply(check_df, is.numeric)
check_df[numeric_col_vec] <- lapply(check_df[numeric_col_vec], signif,
                                    digits = 4)
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(check_df, file.path(csv_dir, "check_permutation_r2.csv"),
                 row.names = FALSE)
print(paste0("Wrote ", file.path(csv_dir, "check_permutation_r2.csv")))
