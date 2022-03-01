rm(list=ls())
library(Seurat)
load("../../../../data/Sydney_stressors/2021-10-28/2021_10_28_Cleaned_Data/final_all_data.RData")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_octoberdata_jackpot1.png"),
                plot1, device = "png", width = 14, height = 8, units = "in")

plot1 <-Seurat::DimPlot(all_data, reduction = "umap",
                        group.by = "OG_condition", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA of October 2021 data"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_octoberdata_umap.png"),
                plot1, device = "png", width = 7, height = 5, units = "in")
