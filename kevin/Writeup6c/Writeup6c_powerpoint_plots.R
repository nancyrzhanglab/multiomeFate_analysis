rm(list=ls())
library(Seurat)
library(Signac)
library(GGally)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat[order(rowSums(tab_mat), decreasing = T)[1:10],]
tab_mat[order(rowSums(tab_mat[,2:4]), decreasing = T)[1:10],]

tab_mat2 <- tab_mat[,1:4]
tab_mat2 <- log10(tab_mat2+1)
tab_mat2 <- apply(as.matrix.noquote(tab_mat2),2,as.numeric)
tab_mat2 <- as.data.frame(tab_mat2)
plot1 <- GGally::ggpairs(tab_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_lineage-count_pairs_log_onlyday10.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

tab_mat2 <- tab_mat[,c(1,5:7)]
tab_mat2 <- log10(tab_mat2+1)
tab_mat2 <- apply(as.matrix.noquote(tab_mat2),2,as.numeric)
tab_mat2 <- as.data.frame(tab_mat2)
plot1 <- GGally::ggpairs(tab_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_lineage-count_pairs_log_onlyweeek5.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

tab_mat2 <- tab_mat
tab_mat2 <- log10(tab_mat2+1)
tab_mat2 <- apply(as.matrix.noquote(tab_mat2),2,as.numeric)
tab_mat2 <- as.data.frame(tab_mat2)
plot1 <- GGally::ggpairs(tab_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_lineage-count_pairs_log.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

###################

tab_mat2 <- tab_mat[,c(1,2,5)]
tab_mat2 <- log10(tab_mat2+1)
tab_mat2 <- apply(as.matrix.noquote(tab_mat2),2,as.numeric)
tab_mat2 <- as.data.frame(tab_mat2)
plot1 <- GGally::ggpairs(tab_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_lineage-count_pairs_log_onlyCIS.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

tab_mat2 <- tab_mat[,c(1,3,6)]
tab_mat2 <- log10(tab_mat2+1)
tab_mat2 <- apply(as.matrix.noquote(tab_mat2),2,as.numeric)
tab_mat2 <- as.data.frame(tab_mat2)
plot1 <- GGally::ggpairs(tab_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_lineage-count_pairs_log_onlyCOCL2.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

tab_mat2 <- tab_mat[,c(1,4,7)]
tab_mat2 <- log10(tab_mat2+1)
tab_mat2 <- apply(as.matrix.noquote(tab_mat2),2,as.numeric)
tab_mat2 <- as.data.frame(tab_mat2)
plot1 <- GGally::ggpairs(tab_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_lineage-count_pairs_log_onlyDABTRAM.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

###################################################

source("../Writeup6b/gene_list.R")

Seurat::Idents(all_data) <- "dataset"
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data2 <- subset(all_data, dataset == "day0")

gene_vec <- sort(c("SOX10", "MITF", "NGFR", 
                   "EGR3", "JUN", "FOSL1", 
                   "CD44", "FN1"))

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
for(treatment in treatment_vec){
  print(treatment)
  
  tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
  lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
  
  ident_vec <- rep(NA, ncol(all_data2))
  names(ident_vec) <- colnames(all_data2)
  cell_names1 <- colnames(all_data2)[which(all_data2$assigned_lineage %in% lineage_names)]
  cell_names2 <- colnames(all_data2)[intersect(which(all_data2$dataset == "day0"),
                                               which(all_data2$assigned_posterior >= 0.5))]
  cell_names_winning <- intersect(cell_names1, cell_names2)
  cell_names_losing <- setdiff(cell_names2, cell_names1)
  ident_vec[cell_names_winning] <- paste0("day0_win_", treatment)
  ident_vec[cell_names_losing] <- paste0("day0_lose_", treatment)
  all_data2$ident <- ident_vec
  Seurat::Idents(all_data2) <- "ident"
  all_data3 <- subset(all_data2, ident %in% c(paste0("day0_win_", treatment), paste0("day0_lose_", treatment)))
  
  pdf(paste0("../../../../out/figures/Writeup6c/Writeup6c_coverage_", treatment, "_5000bp.pdf"),
      onefile = T,
      height = 2.5)
  for(gene in gene_vec){
    plot1 <- Signac::CoveragePlot(
      object = all_data3,
      region = gene,
      features = gene,
      extend.upstream = 5000,
      extend.downstream = 5000
    )
    
    print(plot1)
  }
  
  dev.off() 
}

###############

Seurat::Idents(all_data) <- "dataset"
plot1 <- Seurat::VlnPlot(all_data, features = c("percent.mt"), pt.size = 0)
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_mitochondrial_violin.png"),
                plot1, device = "png", width = 8, height = 4, units = "in")

