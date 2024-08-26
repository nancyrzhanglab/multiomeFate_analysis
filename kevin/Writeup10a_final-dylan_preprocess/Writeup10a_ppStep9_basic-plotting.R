rm(list=ls())
library(Seurat)
library(Signac)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep8_combining.RData"))
treatment_vec <- c("CIS", "COCL2", "DABTRAM")

Seurat::DefaultAssay(all_data) <- "Saver"
plot1 <-Seurat::DimPlot(all_data, 
                        reduction = "saver.umap",
                        group.by = "dataset", 
                        label = TRUE,
                        repel = TRUE, 
                        label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (Saver),\n", length(Seurat::VariableFeatures(all_data[["Saver"]])), " genes for ", length(Seurat::Cells(all_data)), " cells"))
ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_saver_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

for(treatment in treatment_vec){
  plot1 <- Seurat::DimPlot(all_data, 
                          reduction = paste0("ft.", treatment, ".umap"),
                          group.by = "dataset", 
                          label = TRUE,
                          repel = TRUE, 
                          label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (FastTopics): ", treatment, "\n", length(Seurat::Cells(all_data)), " cells"))
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_fasttopic_umap_", treatment, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}