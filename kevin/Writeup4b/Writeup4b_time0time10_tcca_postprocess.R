rm(list=ls())
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_tcca.RData")

library(Seurat); library(Signac)
multiSVD_obj$param
multiSVD_obj$tcca_obj$tilt_perc

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "original_dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Melanoma cells (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_common-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(all_data, reduction = "distinct1_tcca",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Melanoma cells (10x, RNA+ATAC)\nDistinct 1 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_distinct1-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <-Seurat::DimPlot(all_data, reduction = "distinct2_tcca",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Melanoma cells (10x, RNA+ATAC)\nDistinct 2 subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_distinct2-tcca.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")
