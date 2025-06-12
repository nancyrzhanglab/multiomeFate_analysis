# based on https://github.com/nancyrzhanglab/multiomeFate_analysis/blob/emilia/emilia/Final/Fig2/Show_in_vitro_heterogneity.R#L242

rm(list = ls())

library(multiomeFate)
library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"
results_dir <- "~/project/Multiome_fate/out/emilia/task2_correlate_fate_potential_and_features_V2/"

all_data <- multiomeFate:::data_loader(which_files = c("fasttopics"))
scores <- read.csv(paste0(results_dir, 'MITF_AXL_UCell_scores.csv'))
rownames(scores) <- scores$X; scores <- scores[,-1]
stopifnot(all(rownames(scores) == Seurat::Cells(all_data)))

all_data$MITF_Program <- scores[,"MITF_Program"]
all_data$AXL_Program <- scores[,"AXL_Program"]

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
for(treatment in treatment_vec){
  all_data_subset <- subset(all_data, dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
  
  tmp <- all_data_subset$AXL_Program
  tmp[which(all_data_subset$dataset %in% c("day0", paste0("week5_", treatment)))] <- NA
  tmp <- scale(tmp)
  tmp <- pmax(pmin(tmp, 2), -2)
  
  all_data_subset$AXL_Program <- tmp
  
  plot1 <- scCustomize::FeaturePlot_scCustom(all_data_subset, 
                                             features = "AXL_Program",
                                             colors_use = c("#604CC3", "bisque", "#FFA500"),
                                             na_color = "gray90",
                                             na_cutoff = NA,
                                             reduction = paste0("ft.", treatment, ".umap"))
  ggplot2::ggsave(plot1, filename = paste0(plot_folder, "Writeup18_axl_", treatment, ".png"),
                  height = 8, width = 8)
  
  # create an ordering of cells
  # https://github.com/satijalab/seurat/issues/5762
  cell_order_idx <- c(which(is.na(all_data_subset$AXL_Program)),
                      sample(which(!is.na(all_data_subset$AXL_Program))))
  cell_order_names <- Seurat::Cells(all_data_subset)[cell_order_idx]
  plot1 <- scCustomize::FeaturePlot_scCustom(all_data_subset, 
                                             features = "AXL_Program",
                                             colors_use = c("#604CC3", "bisque", "#FFA500"),
                                             na_color = "gray90",
                                             na_cutoff = NA,
                                             order = FALSE,
                                             reduction = paste0("ft.", treatment, ".umap"),
                                             cells = cell_order_names,
                                             pt.size = 2)
  plot1 <- plot1 + Seurat::NoAxes() + Seurat::NoLegend() + ggplot2::ggtitle("")
  ggplot2::ggsave(plot1, filename = paste0(plot_folder, "Writeup18_axl_", treatment, "_cleaned.png"),
                  height = 3, width = 3)
}

for(treatment in treatment_vec){
  all_data_subset <- subset(all_data, dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
  
  plot1 <- scCustomize::DimPlot_scCustom(all_data_subset, 
                                         group.by = "dataset",
                                         reduction = paste0("ft.", treatment, ".umap"), 
                                         colors_use = all_data@misc[["dataset_colors"]])
  ggplot2::ggsave(plot1, filename = paste0(plot_folder, "Writeup18_fasttopic-umap_", treatment, ".png"),
                  height = 8, width = 8)
  
  plot1 <- scCustomize::DimPlot_scCustom(all_data_subset, 
                                         group.by = "dataset",
                                         reduction = paste0("ft.", treatment, ".umap"), 
                                         colors_use = all_data@misc[["dataset_colors"]],
                                         pt.size = 2)
  plot1 <- plot1 + Seurat::NoAxes() + Seurat::NoLegend() + ggplot2::ggtitle("")
  ggplot2::ggsave(plot1, filename = paste0(plot_folder, "Writeup18_fasttopic-umap_", treatment, "_cleaned.png"),
                  height = 3, width = 3)
}
