rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")
source("gene_list.R")

Seurat::DefaultAssay(all_data) <- "chromvar"

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
gene_vec <- sort(unique(unlist(keygenes)))
gene_vec <- gene_vec[which(gene_vec %in% rownames(all_data[["RNA"]]))]

for(treatment in treatment_vec){
  print("====")
  print(treatment)
  
  ident_vec <- all_data$dataset
  names(ident_vec) <- colnames(all_data)
  lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
  cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
  cell_names2 <- colnames(all_data)[which(all_data$dataset == "day0")]
  cell_names_winning <- intersect(cell_names1, cell_names2)
  cell_names_losing <- setdiff(cell_names2, cell_names1)
  ident_vec[cell_names_winning] <- paste0("day0_win_", treatment)
  ident_vec[cell_names_losing] <- paste0("day0_lose_", treatment)
  all_data$ident <- ident_vec
  Seurat::Idents(all_data) <- "ident"
  
  de_res <- Seurat::FindMarkers(
    object = all_data,
    ident.1 = paste0("day0_win_", treatment),
    ident.2 = paste0("day0_lose_", treatment),
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    verbose = F
  )
  
  idx <- which(de_res[,"p_val_adj"] <= 1e-2)
  if(length(idx) == 0) idx <- 1
  print(length(idx))
  
  ## hacked from https://github.com/stuart-lab/signac/blob/master/R/visualization.R
  # plot1 <- Signac::MotifPlot(
  #   object = all_data,
  #   motifs = rownames(de_res)[idx],
  #   assay = "ATAC"
  # )
  motifs <- rownames(de_res)[idx]
  data.use <- Signac::GetMotifData(object = all_data, 
                                    assay = "ATAC", 
                                    slot = "pwm")
  data.use <- data.use[motifs]
  names(data.use) <- Signac::GetMotifData(
    object = all_data, 
    assay = "ATAC", 
    slot = "motif.names"
  )[motifs]
  names(data.use) <- sapply(1:length(names(data.use)), function(i){
    paste0(names(data.use)[i], ", ", rownames(de_res)[idx][i], ", padj=", formatC(de_res[idx[i],"p_val_adj"], format = "e", digits = 2))
  })
  plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
  plot1 <- plot1 + ggplot2::theme_bw()
  
  if(treatment == "DABTRAM"){
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6b/Writeup6b_motif_", treatment, ".png"),
                    plot1, device = "png", width = 15, height = 10, units = "in")
  } else {
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6b/Writeup6b_motif_", treatment, ".png"),
                    plot1, device = "png", width = 10, height = 10, units = "in")
  }
  
}

