rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("gene_list.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::Idents(all_data) <- "dataset"
Seurat::DefaultAssay(all_data) <- "ATAC"

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
  
  pdf(paste0("../../../../out/figures/Writeup6b/Writeup6b_coverage_", treatment, ".pdf"),
      width = 6, height = 5)
  for(gene in gene_vec){
    print(gene)
    
    plot1 <- Signac::CoveragePlot(
      object = all_data,
      region = gene,
      features = gene,
      extend.upstream = 5000,
      extend.downstream = 5000
    )
    plot1
  }
  graphics.off()
}


