rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver.RData")

library(Seurat)
library(Signac)
library(SAVER)

Seurat::DefaultAssay(all_data) <- "Saver"
plot1 <-Seurat::DimPlot(all_data, reduction = "saverumap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (SAVER)"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_rna-saver_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "saverumap",
                             slot = "scale.data",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E", 
                                               "CD44", "LOXL2", "ID3")),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_rna-saver_umap_jackpot1.png"),
                plot1, device = "png", width = 12, height = 12, units = "in")

###############################

