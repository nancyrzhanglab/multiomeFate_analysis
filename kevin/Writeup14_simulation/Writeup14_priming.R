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

coefficient_intercept <- -1.5
coefficient_vec <- rep(0, ncol(cell_features))
names(coefficient_vec) <- colnames(cell_features)
coefficient_vec[1:5] <- seq(0.5, 0.1, length.out = 5)

.compute_mean_total_cells(cell_features, coefficient_intercept, coefficient_vec)

#########

num_lineages <- 50
lineage_spread <- 0.1
lineage_prior <- rep(1/num_lineages, length = num_lineages)

set.seed(10)
simulation_res <- multiomeFate:::generate_simulation(
  embedding_mat = cell_features,
  bool_add_randomness = TRUE,
  coefficient_intercept = coefficient_intercept,
  embedding_coefficient_vec = coefficient_vec,
  lineage_spread = lineage_spread,
  lineage_prior = lineage_prior,
  num_lineages = num_lineages,
  verbose = 3
)

filename <- paste0(plot_folder, "lineage-mean-variance.png")
.plot_mean_variance(filename = filename,
                    simulation_res = simulation_res)



