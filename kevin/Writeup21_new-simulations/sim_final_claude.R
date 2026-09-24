# sim_final_claude.R
# Kevin Z. Lin (drafted by Claude), 2026-09-23
#
# Final-round driver for both Figure 4 sweeps: gene recovery by CYFER, CoSPAR
# and lineage-DE along the t2 Gini axis (Figure 4A) and the within-clone AED
# axis (Figure 4B), at the design settled by the `_r2` trailblazing round
# (simulation_design_claude.md, Section 8). Nothing about the simulation
# changes from `_r2` except the replicate count, which rises from 2 to 20.
#
# Run:  Rscript sim_final_claude.R gini   # Figure 4A
#       Rscript sim_final_claude.R aed    # Figure 4B
# Both write outputs suffixed `_final`, beside the `_r2` trailblazing files.
# The two axes share nothing at run time, so launch them in two shells to run
# them concurrently; each is single-core and takes about 2.2 hours (140
# datasets at roughly 0.9 minutes each).
#
# Dry:   WRITEUP21_DRY=1 Rscript sim_final_claude.R gini   (2 levels x 2
#        replicates, outputs suffixed `_final_dry`; about 4 minutes)
# PCs:   WRITEUP21_PCA=<d> Rscript sim_final_claude.R gini (a non-default
#        embedding dimension; outputs suffixed `_final_pca<d>`)
# Force: WRITEUP21_FORCE=1 allows overwriting an existing RDS of this run.
#
# Watch it with `tail -f` on the progress file, whose path the first lines of
# the run print: every stage of every dataset appends a time-stamped line with
# the elapsed time and an estimate of the time remaining, and each finished
# level appends its mean Spearman per method.

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

# Which axis ------------------------------------------------------------------

arg_vec <- commandArgs(trailingOnly = TRUE)
axis <- if(length(arg_vec) >= 1) arg_vec[1] else ""
if(!axis %in% c("gini", "aed")){
  stop("usage: Rscript sim_final_claude.R <gini|aed>. The axis argument is ",
       "required: `gini` is Figure 4A (levels are t2 Gini targets, the mean ",
       "squared AED held) and `aed` is Figure 4B (levels are mean squared ",
       "AED targets, the t2 Gini held)")
}

# Settings ---------------------------------------------------------------------

bool_dry <- Sys.getenv("WRITEUP21_DRY") == "1"
bool_force <- Sys.getenv("WRITEUP21_FORCE") == "1"
d_pca <- as.numeric(Sys.getenv("WRITEUP21_PCA", unset = "20"))
num_replicates <- 20

# `_final` marks the final round; the trailblazing files ("" = 10 PCs,
# "_pca20", "_r2") stay on disk as the record.
suffix <- paste0("_final",
                 if(d_pca != 20) paste0("_pca", d_pca) else "",
                 if(bool_dry) "_dry" else "")

# Both ladders, both held values and every generator knob are the settled
# `_r2` values (memo Section 8). They live here rather than in two drivers so
# that the two axes cannot drift apart between rounds.
if(axis == "gini"){
  # Gini 0.2 is unreachable at the fixed AED of 0.4 (count noise dominates at
  # its tiny latent scale), so the ladder starts at 0.3.
  level_vec <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)   # t2 Gini targets
  fixed_val <- 0.4                                    # mean squared AED held
  seed_base <- 1000
} else {
  level_vec <- c(0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0)  # mean sq AED targets
  # Memo Section 8: 0.5 is the held value of both trailblazing rounds and
  # stands unless this round's 20 replicates move it. Under the extinct-clone
  # exclusion CoSPAR's `_r2` curve reads its last easy level as 0.6, but on
  # two replicates at levels whose SD reaches 0.3.
  fixed_val <- 0.5                                      # t2 Gini held
  seed_base <- 2000
}
if(bool_dry){
  level_vec <- level_vec[c(1, length(level_vec))]
  num_replicates <- 2
}

rds_file <- file.path(out_dir, paste0("sim_", axis, "_claude", suffix, ".rds"))
progress_file <- file.path(out_dir,
                           paste0("progress_", axis, "_claude", suffix, ".txt"))

# A finished sweep is over two hours of compute, so refuse to start on top of
# one by accident.
if(file.exists(rds_file) && !bool_force){
  stop("`", rds_file, "` already exists. Move it aside, or re-run with ",
       "WRITEUP21_FORCE=1 to overwrite it")
}

# Run --------------------------------------------------------------------------

print(paste0("multiomeFate ", utils::packageVersion("multiomeFate"), ", ",
             R.version.string))

# The seed scheme is unchanged from `_r2` (replicate r of level j uses
# seed_base + 100 j + r, and calibration draws sit at + 50 + k), so replicates
# 1 and 2 regenerate the two trailblazing datasets of each level exactly: a
# free check that nothing in the pipeline drifted, with replicates 3 to 20 new.
result_list <- run_sweep(
  axis = axis,
  level_vec = level_vec,
  fixed_val = fixed_val,
  num_replicates = num_replicates,
  out_dir = out_dir,
  python_path = python_path,
  script_path = file.path(script_dir, "run_cospar.py"),
  rds_file = rds_file,
  progress_file = progress_file,
  d_pca = d_pca,
  kappa = 1.5,
  run_label = paste0(axis, suffix),
  spread_variation = "noncausal",
  tau_delta = 0.6,
  seed_base = seed_base,
  verbose = 1)

print(result_list$summary[, c("level", "method", "spearman_mean",
                              "spearman_sd", "pearson_mean", "jaccard_mean",
                              "gini_t2_mean", "aed_mean_mean")])
print(paste0("wrote ", rds_file))
