rm(list=ls())
# meant to be run from the laptop currently
library(Seurat)
source("../multiomeFate_analysis/kevin/Writeup3d/funcs.R")
source("../multiomeFate_analysis/kevin/Writeup3e/select_cells.R")

load("../../dbox_MultiomeFate/data/ShafferLab/DLS_scRNA-seq_for_Zhang_Lab/all_data_SCT.RData")
