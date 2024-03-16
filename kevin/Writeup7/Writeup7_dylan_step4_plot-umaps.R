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

########

tab_mat <- table(all_data$assigned_lineage, all_data$OG_condition)
tab_mat <- log10(tab_mat+1)

df <- as.data.frame(apply(as.matrix.noquote(tab_mat),2,as.numeric))
p1 <- GGally::ggpairs(df,
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup7/Writeup7_dylan_step4_lineage-table.png",
                p1, device = "png", width = 8, height = 8, units = "in")



