rm(list=ls())
load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed_de.RData")

library(Seurat); library(Signac)

peak_names <- sort(unique(mbrain3[["ATAC"]]@links$peak))
gene_names <- sort(unique(mbrain3[["ATAC"]]@links$gene))

mbrain3[["pca"]] <- NULL
mbrain3[["lsi"]] <- NULL
mbrain3[["umap.rna"]] <- NULL
mbrain3[["umap.atac"]] <- NULL
mbrain3[["wnn.umap"]] <- NULL

Seurat::DefaultAssay(mbrain3) <- "SCT"
mbrain3 <- Seurat::RunPCA(mbrain3, features = gene_names, verbose = F)
set.seed(10)
mbrain3 <- Seurat::RunUMAP(mbrain3, dims = 1:50, reduction.name = 'umap.rna', 
                           reduction.key = 'rnaUMAP_')

#######

Seurat::DefaultAssay(mbrain3) <- "ATAC"
mbrain3[["ATAC"]]@var.features <- peak_names
mbrain3 <- Signac::RunSVD(mbrain3, features = mbrain3[["ATAC"]]@var.features)
set.seed(10)
mbrain3 <- Seurat::RunUMAP(mbrain3, reduction = 'lsi', dims = 2:50, 
                           reduction.name = "umap.atac", reduction.key = "atacUMAP_")

set.seed(10)
mbrain3 <- Seurat::FindMultiModalNeighbors(mbrain3, reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:50, 2:50))
mbrain3 <- Seurat::RunUMAP(mbrain3, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                reduction.key = "wnnUMAP_")

##############

plot1 <- Seurat::DimPlot(mbrain3, reduction = "umap.rna", group.by = "new_seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nRNA (Seurat clusters)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de_rna_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain3, reduction = "umap.rna", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nRNA (Celltype)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de_rna.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

##

plot1 <- Seurat::DimPlot(mbrain3, reduction = "umap.atac", group.by = "new_seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nATAC (Seurat clusters)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de_atac_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain3, reduction = "umap.atac", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nATAC (Celltype)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de_atac.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

##

plot1 <- Seurat::DimPlot(mbrain3, reduction = "wnn.umap", group.by = "new_seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nWNN (Seurat clusters)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de_wnn_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(mbrain3, reduction = "wnn.umap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nWNN (Celltype)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de_wnn.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")
