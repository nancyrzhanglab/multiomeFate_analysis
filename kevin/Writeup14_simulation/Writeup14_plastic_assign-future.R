rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_plastic-setting_"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_plastic-setting_"
func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"

source(paste0(func_folder, "func.R"))
source(paste0(func_folder, "func_priming.R"))
source(paste0(func_folder, "func_assign-future-cells.R"))
source(paste0(func_folder, "func_seurat.R"))

#######

all_data <- multiomeFate:::data_loader(which_files = c( "rna", "saver", "fasttopics"))
keep_vec <- rep(FALSE, length(Seurat::Cells(all_data)))
keep_vec[all_data$dataset %in% c("day10_COCL2", "week5_COCL2")] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)
all_data[["fasttopic.CIS"]] <- NULL
all_data[["ft.CIS.umap"]] <- NULL
all_data[["fasttopic.DABTRAM"]] <- NULL
all_data[["ft.DABTRAM.umap"]] <- NULL

load(paste0(out_folder, "simulation.RData"))

# count how many future cells there's supposed to be
total_future_cells <- sum(simulation_res$lineage_future_size)
current_future_cells <- length(which(all_data$dataset == "week5_COCL2"))

total_future_cells; current_future_cells

# here, since total_future_cells < current_future_cells, we simply throw out cells
simulation_data <- .form_simulation_seurat_fate(final_fit = final_fit,
                                                simulation_res = simulation_res)

# subset only "day10_COCL2" and "week5_COCL2"
lineage_vec <- all_data$assigned_lineage
lineage_vec[Seurat::Cells(simulation_data)] <- simulation_data$assigned_lineage
simulation_data$assigned_lineage <- lineage_vec

set.seed(10)
all_data <- .assign_future_cells(
  seurat_obj = all_data,
  future_ident = "week5_COCL2",
  ident = "dataset",
  lineage_future_size = simulation_res$lineage_future_size,
  lineage_variable = "assigned_lineage"
)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
sum(abs(sort(tab_mat[,2])-sort(simulation_res$lineage_future_size)))


