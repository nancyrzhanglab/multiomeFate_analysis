rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6b/Writeup6b_DE_day10_expanding.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::Idents(all_data) <- "dataset"
Seurat::DefaultAssay(all_data) <- "ATAC"

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

for(kk in 2:length(treatment_vec)){
  treatment <- treatment_vec[kk]
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
  
  gene_vec1 <- rownames(de_list[[kk]])[which(de_list[[kk]][,"p_val_adj"] <= 1e-2)]
  gene_vec2 <- rownames(de_list[[kk]])[order(de_list[[kk]][,"p_val_adj"], decreasing = F)[1:50]]
  gene_vec <- intersect(gene_vec1, gene_vec2)
  tmp <- grep("-", gene_vec)
  if(length(tmp) > 0) gene_vec <- gene_vec[-tmp]
  tmp <- sapply(gene_vec, function(gene){
    zz <- tryCatch({Signac:::FindRegion(
      object = all_data,
      region = gene,
      sep = c("-", "-"),
      extend.upstream = 1000,
      extend.downstream = 1000
    )}, error = function(e){NA})
    if(!is.na(zz)) TRUE else FALSE
  })
  gene_vec <- gene_vec[tmp]
  print(length(gene_vec))
  
  if(length(gene_vec) > 0){
    pdf(paste0("../../../../out/figures/Writeup6b/Writeup6b_coverage_", treatment, "_DEgenes_1000bp.pdf"),
        onefile = T)
    for(gene in gene_vec){
      plot1 <- Signac::CoveragePlot(
        object = all_data,
        region = gene,
        features = gene,
        extend.upstream = 1000,
        extend.downstream = 1000
      )
      
      print(plot1)
    }
    
    dev.off() 
  }
}



