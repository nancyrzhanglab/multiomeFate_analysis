rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[tab_mat[,"day0"] > 0,]
lineage_names <- 

# for each lineage, compute