rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../data/Greenleaf_humancortex/09282021_seurat_processed.RData")

plot1 <- Seurat::DimPlot(cortex, reduction = "umap.rna",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human cortex (RNA)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_greenleaf_rna.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(cortex, reduction = "umap.atac",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human cortex (ATAC)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_greenleaf_atac.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(cortex, reduction = "wnn.umap",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("Human cortex (WNN)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_greenleaf_wnn.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")