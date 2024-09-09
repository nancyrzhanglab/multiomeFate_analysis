rm(list=ls())
library(multiomeFate)
library(Seurat)
# library(ggpubr)

source('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/git/multiomeFate/R/data_loader.R')
# =============================================================================
# Load data
# =============================================================================
# args = commandArgs(trailingOnly=TRUE)
# TIME = args[1]
# TREATMENT = args[2]
TIME = 'day10' # 'week5', 'day10', or 'day0'
TREATMENT = 'DABTRAM' # 'COCL2', 'DABTRAM', or 'CIS'
SAMPLE_NAME = paste0(TIME, '_', TREATMENT)

# all_data = multiomeFate:::data_loader(which_files = c("saver"))
load('~/Downloads/Writeup10a_data_empty.RData')
load('~/Downloads/Writeup10a_data_saver.RData')
all_data[['Saver']] <- all_data_saver
output_dir = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/processed_data/"
output_dir = "/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/"
# =============================================================================
# Wrangle data
# =============================================================================
all_data_use = subset(all_data, dataset == SAMPLE_NAME)

Seurat::DefaultAssay(all_data_use) <- "Saver"
all_data_use <- NormalizeData(all_data_use)
all_data_use <- FindVariableFeatures(all_data_use)
all_data_use <- ScaleData(all_data_use)
all_data_use <- RunPCA(all_data_use)
all_data_use <- RunUMAP(all_data_use, dims = 1:30)

data_pca <- all_data_use@reductions[["pca"]]@cell.embeddings

write.csv(data_pca, paste0(output_dir, "data_pca_saver_sample_", SAMPLE_NAME, ".csv"))

