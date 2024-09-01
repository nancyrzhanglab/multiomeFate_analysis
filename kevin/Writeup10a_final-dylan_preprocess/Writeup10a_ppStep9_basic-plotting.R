rm(list=ls())
library(Seurat)
library(Signac)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"

#######

all_data <- data_loader(which_files = c("rna", "saver", "peakvi", "fasttopics", "wnn"))

dimred_vec <- names(all_data@reductions)
dimred_vec <- dimred_vec[grep("umap", dimred_vec)]

pdf(paste0(plot_folder, "Writeup10a_all-umaps.pdf"),
    onefile = TRUE, width = 8, height = 5)

for(kk in 1:length(dimred_vec)){
  plot1 <- scCustomize::DimPlot_scCustom(all_data, 
                           reduction = dimred_vec[kk],
                           group.by = "dataset",
                           colors_use = all_data@misc$dataset_colors)
  plot1 <- plot1 + ggplot2::ggtitle(dimred_vec[kk])
  print(plot1)
}

dev.off()