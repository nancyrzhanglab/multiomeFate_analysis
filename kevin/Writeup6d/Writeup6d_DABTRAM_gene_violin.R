rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################################

treatment <- "DABTRAM"
Seurat::DefaultAssay(all_data) <- "RNA"

treatment <- "DABTRAM"
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
surviving_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
dying_lineages <- rownames(tab_mat)[which(apply(tab_mat,1,max)<=1)]
winning_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% surviving_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
dying_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% dying_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
length(winning_idx); length(dying_idx)
ident_vec <- rep(NA, ncol(all_data))
ident_vec[winning_idx] <- paste0("day0_win_", treatment)
ident_vec[dying_idx] <- paste0("day0_lose_", treatment)
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"
table(Seurat::Idents(all_data))

all_data <- subset(all_data, ident %in% c(paste0("day0_win_", treatment), paste0("day0_lose_", treatment)))

######################

source("../Writeup6b/gene_list.R")
source("gene_list_csc.R")
important_genes <- sort(unique(c(unlist(keygenes), keygenes_csc)))

gene_vec <- sort(unique(c(important_genes, all_data[["Saver"]]@var.features)))
gene_vec <- intersect(gene_vec, rownames(all_data[["RNA"]]@counts))

important_genes <- intersect(important_genes, gene_vec)

######################

set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  features = gene_vec,
  ident.1 = paste0("day0_win_", treatment),
  ident.2 = paste0("day0_lose_", treatment),
  logfc.threshold = 0,
  min.pct = 0
)
de_res[c("JUN","JUNB","FOSL1"),]

for(k in 1:3){
  print(k)
  genes <- important_genes[((k-1)*25+1):min(k*25, length(important_genes))]
  
  Seurat::DefaultAssay(all_data) <- "RNA"
  # see https://github.com/satijalab/seurat/issues/312
  plot1 <- Seurat::VlnPlot(all_data, 
                           features = genes,
                           slot = "data",
                           ncol = 5,
                           pt.size = 0.5, combine = F)
  title_vec <- sapply(genes, function(x){
    paste0(x, ": -Log10 Pvalue: ", round(-log10(de_res[x,"p_val"]),2), "\nDiff: ", round(de_res[x,"avg_log2FC"],2))
  })
  for(i in 1:length(title_vec)){
    plot1[[i]] <- plot1[[i]] + ggplot2::ggtitle(title_vec[i]) + ggplot2::theme_bw()
    plot1[[i]] <- plot1[[i]] + Seurat::NoLegend()
  }
 
  plot_all <- cowplot::plot_grid(plotlist = plot1, nrow = 5)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6d/Writeup6d_DABTRAM_genes",k, "_RNA-data_violinplots.png"),
                  plot_all, device = "png", width = 20, height = 20, units = "in")
}

#################

Seurat::DefaultAssay(all_data) <- "Saver"
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  features = intersect(gene_vec, all_data[["Saver"]]@var.features),
  ident.1 = paste0("day0_win_", treatment),
  ident.2 = paste0("day0_lose_", treatment),
  logfc.threshold = 0,
  min.pct = 0
)
de_res[c("JUN","JUNB","FOSL1"),]

important_genes2 <- sort(intersect(important_genes, all_data[["Saver"]]@var.features))
for(k in 1:2){
  print(k)
  genes <- important_genes2[((k-1)*25+1):min(k*25, length(important_genes2))]
  
  Seurat::DefaultAssay(all_data) <- "Saver"
  plot1 <- Seurat::VlnPlot(all_data, 
                           features = genes,
                           slot = "data",
                           ncol = 5,
                           pt.size = 0.5, combine = F)
  title_vec <- sapply(genes, function(x){
    paste0(x, ": -Log10 Pvalue: ", round(-log10(de_res[x,"p_val"]),2), "\nDiff: ", round(de_res[x,"avg_log2FC"],2))
  })
  for(i in 1:length(title_vec)){
    plot1[[i]] <- plot1[[i]] + ggplot2::ggtitle(title_vec[i]) + ggplot2::theme_bw()
    plot1[[i]] <- plot1[[i]] + Seurat::NoLegend()
  }
  
  plot_all <- cowplot::plot_grid(plotlist = plot1, nrow = 5)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6d/Writeup6d_DABTRAM_genes",k, "_RNA-Saver_violinplots.png"),
                  plot_all, device = "png", width = 20, height = 20, units = "in")
}
