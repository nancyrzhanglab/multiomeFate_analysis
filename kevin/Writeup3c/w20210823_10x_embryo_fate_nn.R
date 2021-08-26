rm(list=ls())
load("../../../../out/kevin/Writeup3c/20210823_10x_embryo_result.RData")

library(Seurat); library(Signac); library(igraph)
source("nn_plot.R")

end_states <- list(16, 6, c(1,2,4), 9)
names(end_states) <- c("Oligo.", "Forebrain", "Cortical1", "Cortical2")
celltype <- mbrain3@meta.data$new_seurat_clusters

#################

tmp <- nn_plot(mbrain3, res$nn_mat, celltype, end_states)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_nndist.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

#############

tmp <- dirnn_plot(mbrain3, res$list_diagnos, celltype, end_states)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_dirnndist.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")




