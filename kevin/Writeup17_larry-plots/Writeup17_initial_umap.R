rm(list=ls())
library(Seurat)
library(ggplot2)

load("/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup13/Writeup13_larry-dataset.RData")
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17/"

col_vec <- c(
  "Monocyte-4" = "#A6CEE3",
  # "Monocyte-4"= "#1F78B4",
  "Monocyte-6" = "#4472C4",
  "Neutrophil-4"= "#FCA5A5",
  # "Neutrophil-4" = "#E31A1C",
  "Neutrophil-6" = "#D55E00",
  # "Undifferentiated-4" = "#D3D3D3",
  "Undifferentiated-4" = "#A9A9A9",
  "Undifferentiated-6" = "#7F7F7F"
)

seurat_object_subset <- subset(seurat_object, Time.point %in% c("4","6"))

plot1 <- scCustomize::DimPlot_scCustom(seurat_object_subset, 
                              reduction = "SPRING", 
                              group.by = "time_celltype",
                              colors_use = col_vec)
ggplot2::ggsave(plot1, filename = paste0(plot_folder, "initial_umap.png"),
                height = 8, width = 8)

# make cleaned version
plot1 <- plot1 + theme_void() +  # Removes all background, axes, and grid lines
  theme(
    plot.title = element_blank(),  # Remove title
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks = element_blank()  # Remove axis ticks
  ) + Seurat::NoLegend() 

# Display the plot
ggplot2::ggsave(plot1, filename = paste0(plot_folder, "initial_umap_cleaned.png"),
                height = 2, width = 2)

