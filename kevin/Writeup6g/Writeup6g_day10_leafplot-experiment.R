rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

gene_vec <- intersect(gene_vec, rownames(all_data[["geneActivity"]]@data))
gene_vec <- intersect(gene_vec, rownames(all_data[["Saver"]]@data))
gene_vec <- sort(gene_vec)
length(gene_vec)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
treatment_vec <- c("CIS", "COCL2", "DABTRAM")

for(treatment in treatment_vec){
  print("====")
  print(treatment)
  
  # find the winning and losing cells
  tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
  tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 5),
                                                which(tab_mat[,paste0("week5_", treatment)] <= 25))]
  tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
  length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)
  
  tier1_idx <- intersect(
    intersect(which(all_data$assigned_lineage %in% tier1_lineages),
              which(all_data$assigned_posterior >= 0.5)),
    which(all_data$dataset == paste0("day10_", treatment))
  )
  tier2_idx <- intersect(
    intersect(which(all_data$assigned_lineage %in% tier2_lineages),
              which(all_data$assigned_posterior >= 0.5)),
    which(all_data$dataset == paste0("day10_", treatment))
  )
  tier3_idx <- intersect(
    intersect(which(all_data$assigned_lineage %in% tier3_lineages),
              which(all_data$assigned_posterior >= 0.5)),
    which(all_data$dataset == paste0("day10_", treatment))
  )
  
  keep_vec <- rep(NA, ncol(all_data))
  keep_vec[tier1_idx] <- paste0("3high_winner_", treatment)
  keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
  keep_vec[tier3_idx] <- paste0("1loser_", treatment)
  table(keep_vec)
  all_data2 <- all_data
  all_data2$keep <- keep_vec
  all_data2 <- subset(all_data2, keep %in% c(paste0("3high_winner_", treatment),
                                             paste0("2mid_winner_", treatment),
                                             paste0("1loser_", treatment)))
  Seurat::DefaultAssay(all_data2) <- "ATAC"
  Seurat::Idents(all_data2) <- "keep"
  
  color_vec <- rep(NA, ncol(all_data2))
  color_palette <- grDevices::colorRampPalette(c("lightgray", "dodgerblue4"))(3)
  color_vec[which(all_data2$keep ==  paste0("1loser_", treatment))] <- color_palette[1]
  color_vec[which(all_data2$keep ==  paste0("2mid_winner_", treatment))] <- color_palette[2]
  color_vec[which(all_data2$keep ==  paste0("3high_winner_", treatment))] <- color_palette[3]
  
  png("../../../../out/figures/Writeup6f/")
  
  
  dev.off() 
}