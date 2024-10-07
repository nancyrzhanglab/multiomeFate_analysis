rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_plastic-setting_"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_priming-setting_"

#######

all_data <- multiomeFate:::data_loader(which_files = "fasttopics")

set.seed(10)
rna_mat <- all_data[["fasttopic.COCL2"]]@cell.embeddings
cell_features <- .preprocess_rna(rna_mat, "day10_COCL2")

set.seed(10)
tmp <- .search_for_priming_parameters(cell_features,
                                      min_cells = 6000,
                                      min_maximum = 10,
                                      verbose = 1)
coefficient_intercept <- tmp$coefficient_intercept
coefficient_vec <- tmp$coefficient_vec
num_cells <- .compute_mean_total_cells(cell_features, coefficient_intercept, coefficient_vec, return_sum = FALSE)
round(quantile(num_cells), 3)

#########

num_lineages <- 50

set.seed(10)
simulation_res <- .assign_lineages_priming(
  embedding_mat = cell_features,
  coefficient_intercept = coefficient_intercept,
  coefficient_vec = coefficient_vec,
  num_lineages = num_lineages,
  num_rounds = 50,
  verbose = 1
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




