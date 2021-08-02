rm(list=ls())
load("../../../../out/kevin/Writeup3b/mouseicb_fate_prep.RData")

# create a new Seurat object
myeloid2 <- Seurat::CreateSeuratObject(counts = Matrix::t(mat_y))
tmp <- Seurat::CreateAssayObject(counts = Matrix::t(mat_x))
myeloid2[["ATAC"]] <- tmp

Seurat::DefaultAssay(myeloid2) <- "RNA"
set.seed(10)
myeloid2 <- Seurat::SCTransform(myeloid2)
myeloid2 <- Seurat::RunPCA(myeloid2, verbose = FALSE)
set.seed(10)
myeloid2 <- Seurat::RunUMAP(myeloid2, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

Seurat::DefaultAssay(myeloid2) <- "ATAC"
set.seed(10)
myeloid2 <- Signac::RunTFIDF(myeloid2)
myeloid2[["ATAC"]]@var.features <- rownames(myeloid2[["ATAC"]]@counts)
myeloid2 <- Signac::RunSVD(myeloid2)
set.seed(10)
myeloid2 <- Seurat::RunUMAP(myeloid2, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

set.seed(10)
myeloid2 <- Seurat::FindMultiModalNeighbors(myeloid2, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
myeloid2 <- Seurat::RunUMAP(myeloid2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
myeloid2[["celltype"]] <- myeloid@meta.data$celltype[cell_idx]

##################

plot1 <- Seurat::DimPlot(myeloid2, reduction = "umap.rna", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse ICB\nRNA")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_seurat_rna.png", 
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(myeloid2, reduction = "umap.atac", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse ICB\nATAC")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_seurat_atac.png", 
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(myeloid2, reduction = "wnn.umap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse ICB\nWNN")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_seurat_wnn.png", 
                plot1, device = "png", width = 5, height = 5, units = "in")

save(myeloid2, file = "../../../../out/kevin/Writeup3b/20210731_mouseicb_preprocess2.RData")

##########################

cell_idx2 <- which(myeloid2@meta.data$celltype == "B16_dICB")
myeloid3 <- Seurat::CreateSeuratObject(counts = myeloid2[["RNA"]]@counts[,cell_idx2])
Seurat::DefaultAssay(myeloid3) <- "RNA"
set.seed(10)
myeloid3 <- Seurat::SCTransform(myeloid3)
myeloid3 <- Seurat::RunPCA(myeloid3, verbose = FALSE)
myeloid3 <- Seurat::FindNeighbors(myeloid3, dims = 1:10)
myeloid3 <- Seurat::FindClusters(myeloid3, resolution = 0.5)

tmp <- rep("Unknown", ncol(myeloid2))
for(i in 1:nrow(myeloid3@meta.data)){
  tmp[which(rownames(myeloid2@meta.data) == rownames(myeloid3@meta.data)[i])] <- paste0("Cluster", myeloid3@meta.data$seurat_clusters[i])
}
myeloid2[["celltype2"]] <- tmp

plot1 <- Seurat::DimPlot(myeloid2, reduction = "umap.rna", group.by = "celltype2", label = TRUE,
                         repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse ICB\nRNA")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_seurat_rna2.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(myeloid2, reduction = "umap.atac", group.by = "celltype2", label = TRUE,
                         repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse ICB\nATAC")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_seurat_atac2.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(myeloid2, reduction = "wnn.umap", group.by = "celltype2", label = TRUE,
                         repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse ICB\nWNN")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_seurat_wnn2.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

save(myeloid2, file = "../../../../out/kevin/Writeup3b/20210731_mouseicb_preprocess2.RData")
