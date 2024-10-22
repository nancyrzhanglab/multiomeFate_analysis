rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup14/Writeup14_priming-setting_"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup14/Writeup14_priming-setting_"
func_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/kevin/Writeup14_simulation/"

source(paste0(func_folder, "func.R"))
source(paste0(func_folder, "func_priming.R"))

#######

all_data <- multiomeFate:::data_loader(which_files = "fasttopics")
load(paste0(out_folder, "simulation.RData"))

# count how many future cells there's supposed to be
total_future_cells <- sum(simulation_res$lineage_future_size)
current_future_cells <- length(which(all_data$dataset == "week5_COCL2"))

total_future_cells; current_future_cells

# here, since total_future_cells > current_future_cells, we randomly combine cells