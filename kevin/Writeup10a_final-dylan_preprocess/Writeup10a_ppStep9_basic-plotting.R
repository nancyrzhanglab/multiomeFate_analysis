rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"

#######

all_data <- multiomeFate:::data_loader(which_files = c("rna", "saver", "peakvi", "fasttopics", "wnn"))

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

##########

all_data$keep <- as.numeric(is.na(all_data$assigned_lineage))

table(all_data$dataset, is.na(all_data$assigned_lineage))

dimred_vec <- names(all_data@reductions)
dimred_vec <- dimred_vec[grep("umap", dimred_vec)]

pdf(paste0(plot_folder, "Writeup10a_all-umaps_diagnostic-missing-lineage.pdf"),
    onefile = TRUE, width = 8, height = 5)

for(kk in 1:length(dimred_vec)){
  plot1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                             reduction = dimred_vec[kk],
                                             features = "keep")
  plot1 <- plot1 + ggplot2::ggtitle(dimred_vec[kk])
  print(plot1)
}

dev.off()
