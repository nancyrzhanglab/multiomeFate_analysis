rm(list=ls())
library(Seurat)

load("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step3_fasttopics.RData")
set.seed(10)

plot1 <- Seurat::DimPlot(all_data, 
                      group.by = "OG_condition",
                      reduction = "saverumap")
plot1 <- plot1 + ggplot2::ggtitle("Saver embedding")
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup7/Writeup7_dylan_step4_umap-saver.png",
                plot1, device = "png", width = 8, height = 7, units = "in")

plot1 <- Seurat::DimPlot(all_data, 
                         group.by = "OG_condition",
                         reduction = "umap")
plot1 <- plot1 + ggplot2::ggtitle("RNA embedding")
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup7/Writeup7_dylan_step4_umap.png",
                plot1, device = "png", width = 8, height = 7, units = "in")


