rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

library(ggplot2)                     
library(GGally)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[tab_mat[,"day0"] > 0,]
lineage_names <- rownames(tab_mat)

data.use <- Signac::GetMotifData(object = all_data,
                                 assay = "ATAC",
                                 slot = "pwm")
data.use <- data.use[motifs]
names(data.use) <- Signac::GetMotifData(
  object = all_data,
  assay = "ATAC",
  slot = "motif.names"
)[motifs]
name_vec <- names(data.use) 

# for each lineage, compute the mean scores
chromvar_scores