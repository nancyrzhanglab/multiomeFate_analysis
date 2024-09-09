rm(list=ls())
library(multiomeFate)
library(dplyr)
library(ggplot2)
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

all_data = multiomeFate:::data_loader(which_files = c("saver"))

output_dir = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/processed_data/"

# =============================================================================
# Wrangle data
# =============================================================================
all_data = subset(all_data, dataset == SAMPLE_NAME)

Seurat::DefaultAssay(data) <- "Saver"
data <- NormalizeData(data)
data <- FindVariableFeatures(data)
data <- ScaleData(data)
data <- RunPCA(data)
data <- RunUMAP(data, dims = 1:30)

data_pca <- data@reductions[["pca"]]@cell.embeddings

write.csv(data_pca, paste0(output_dir, "data_pca_saver_sample_", SAMPLE_NAME, ".csv"))
