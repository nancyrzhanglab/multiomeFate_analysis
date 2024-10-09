rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_plastic-setting_"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_plastic-setting_"
func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"

source(paste0(func_folder, "func.R"))
source(paste0(func_folder, "func_priming.R"))

#######

all_data <- multiomeFate:::data_loader(which_files = "fasttopics")

set.seed(10)
rna_mat <- all_data[["fasttopic.COCL2"]]@cell.embeddings
cell_features <- .preprocess_rna(rna_mat, "day10_COCL2")

coefficient_intercept <- -1.5
coefficient_vec <- rep(0, ncol(cell_features))
names(coefficient_vec) <- colnames(cell_features)
coefficient_vec[1:5] <- seq(0.5, 0.1, length.out = 5)

.compute_mean_total_cells(cell_features, coefficient_intercept, coefficient_vec, return_sum = TRUE)

#########

num_lineages <- 50

set.seed(10)
simulation_res <- multiomeFate:::generate_simulation_plastic(
  embedding_mat = cell_features,
  bool_add_randomness = TRUE,
  coefficient_intercept = coefficient_intercept,
  embedding_coefficient_vec = coefficient_vec,
  num_lineages = num_lineages,
  lineage_mean_spread = 1, 
  lineage_sd_spread = NA,
  verbose = 3
)

filename <- paste0(plot_folder, "lineage-mean-variance.png")
.plot_mean_variance(filename = filename,
                    simulation_res = simulation_res)

##################

# required ingredients:
# cell_features
# cell_lineage (from simulation_res$lineage_assignment)
# lineage_future_count (from simulation_res$lineage_future_size)
# tab_mat (from simulation_res$lineage_assignment and lineage_future_count)

# try out our fate potential method
cell_lineage <- as.character(simulation_res$lineage_assignment)
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- simulation_res$lineage_future_size
tmp <- table(simulation_res$lineage_assignment)
lineage_current_count <- as.numeric(tmp); names(lineage_current_count) <- names(tmp)
lineage_current_count <- lineage_current_count[names(lineage_future_count)]
tab_mat <- cbind(lineage_current_count, lineage_future_count)
colnames(tab_mat) <- c("now", "future")

#################
# start cross validation

set.seed(10)
fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = "future",
  lineage_future_count = lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 10,
  tab_mat = tab_mat,
  num_folds = 10,
  verbose = 2
)

final_fit <- multiomeFate:::lineage_cv_finalize(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  fit_res = fit_res,
  lineage_future_count = lineage_future_count
)

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(fit_res, final_fit, simulation_res,
     date_of_run, session_info,
     file = paste0(out_folder, "simulation.RData"))

print("Done! :)")
