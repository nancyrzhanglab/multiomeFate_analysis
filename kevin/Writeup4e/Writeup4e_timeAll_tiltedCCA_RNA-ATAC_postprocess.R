rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_tcca_RNA-ATAC.RData")

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_tcca_RNA-ATAC_umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "Phase", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nCommon subspace, Cell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_tcca_RNA-ATAC_umap_common_cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

############

plot2 <- Seurat::DimPlot(all_data, reduction = "distinct1_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nRNA (Saver) distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_tcca_RNA-ATAC_umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(all_data, reduction = "distinct1_tcca",
                         group.by = "Phase", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nRNA distinct subspace, Cell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_tcca_RNA-ATAC_umap_distinct1_cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

############

plot3 <- Seurat::DimPlot(all_data, reduction = "distinct2_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nATAC distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_tcca_RNA-ATAC_umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(all_data, reduction = "distinct2_tcca",
                         group.by = "Phase", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nATAC distinct subspace, Cell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_tcca_RNA-ATAC_umap_distinct2_cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

######################

set.seed(10)
all_data <- Seurat::FindMultiModalNeighbors(
  all_data, reduction.list = list("saverpca", "lsi"), 
  dims.list = list(1:50, 2:50), modality.weight.name = "Saver.weight"
)
all_data <- Seurat::RunUMAP(all_data, 
                            nn.name = "weighted.nn", 
                            reduction.name = "wnn.umap", 
                            reduction.key = "wnnUMAP_")

plot3 <- Seurat::DimPlot(all_data, reduction = "wnn.umap",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("WNN: RNA (Saver)+ATAC"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_wnn_RNA-ATAC_umap.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(all_data, reduction = "wnn.umap",
                         group.by = "Phase", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("WNN: RNA (Saver)+ATAC\nCell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_wnn_RNA-ATAC_umap_cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###################

plot1 <- Seurat::DimPlot(all_data, reduction = "saverumap",
                         group.by = "Phase", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (Saver)\nCell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_saver_cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(all_data, reduction = "atac.umap",
                         group.by = "Phase", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (TF-IDF)\nCell-cycle phase"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_atac_cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


