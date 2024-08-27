rm(list=ls())
library(Seurat)
library(Signac)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc-step2.RData"))

#######

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[order(rowSums(tab_mat),decreasing = TRUE),]
df <- as.data.frame(apply(as.matrix.noquote(log10(tab_mat+1)),2,as.numeric))
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points",
                                                             alpha = 0.2, 
                                                             shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_lineage_pairsplot.png"),
                p1, device = "png", width = 8, height = 8, units = "in")
