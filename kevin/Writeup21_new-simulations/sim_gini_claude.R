# sim_gini_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-22
#
# Figure 4A driver: gene recovery by CYFER, CoSPAR and lineage-DE as the t2
# Gini of clone sizes rises, with the mean squared AED held fixed
# (simulation_design_claude.md, Sections 3, 5 and 8). Trailblazing round:
# seven Gini targets x two replicates at the provisional fixed AED of 0.4
# (h2 = 0.6) and tau_delta = 0.6.
#
# Run:   Rscript sim_gini_claude.R
# Dry:   WRITEUP21_DRY=1 Rscript sim_gini_claude.R   (2 levels x 1 replicate,
#        outputs suffixed `_dry`)
# PCs:   WRITEUP21_PCA=20 Rscript sim_gini_claude.R  (20-PC shared embedding
#        instead of 10; outputs suffixed `_pca20`)
# Progress is appended to OUT_ROOT/Writeup21_new-simulations/
# progress_gini_claude.txt; the RDS lands beside it after every level.

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
python_path <- "/opt/miniconda3/envs/cospar/bin/python"

source(file.path(script_dir, "func_generate_claude.R"))
source(file.path(script_dir, "func_methods_claude.R"))
source(file.path(script_dir, "func_sweep_claude.R"))
source(file.path(script_dir, "cospar_flat_io.R"))

# Settings ---------------------------------------------------------------------

bool_dry <- Sys.getenv("WRITEUP21_DRY") == "1"
d_pca <- as.numeric(Sys.getenv("WRITEUP21_PCA", unset = "10"))
suffix <- paste0(if(d_pca != 10) paste0("_pca", d_pca) else "",
                 if(bool_dry) "_dry" else "")

level_vec <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9)   # t2 Gini targets
fixed_aed <- 0.4                                     # mean squared AED held
num_replicates <- 2
if(bool_dry){
  level_vec <- level_vec[c(1, length(level_vec))]
  num_replicates <- 1
}

# Run --------------------------------------------------------------------------

result_list <- run_sweep(
  axis = "gini",
  level_vec = level_vec,
  fixed_val = fixed_aed,
  num_replicates = num_replicates,
  out_dir = out_dir,
  python_path = python_path,
  script_path = file.path(script_dir, "run_cospar.py"),
  rds_file = file.path(out_dir, paste0("sim_gini_claude", suffix, ".rds")),
  progress_file = file.path(out_dir,
                            paste0("progress_gini_claude", suffix, ".txt")),
  d_pca = d_pca,
  kappa = 1.5,
  run_label = paste0("gini", suffix),
  spread_variation = "noncausal",
  tau_delta = 0.6,
  seed_base = 1000,
  verbose = 1)

print(result_list$summary[, c("level", "method", "spearman_mean",
                              "pearson_mean", "jaccard_mean",
                              "gini_t2_mean", "aed_mean_mean")])
