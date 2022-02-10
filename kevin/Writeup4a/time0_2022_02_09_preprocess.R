rm(list=ls())
library(Seurat)
load("../../../../out/kevin/Writeup4a/2022-02-09_time0_formatted.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

plot1 <- Seurat::VlnPlot(t0_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_QC_metric.png"),
                plot1, device = "png", width = 10, height = 6, units = "in")

plot1 <- Seurat::FeatureScatter(t0_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- plot1 + Seurat::NoLegend()
plot2 <- Seurat::FeatureScatter(t0_obj, feature1 = "nCount_Lineage", feature2 = "nFeature_Lineage")
plot2 <- plot2 + Seurat::NoLegend()
plot1 <- plot1 + plot2
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_QC_scatterplot.png"),
                plot1, device = "png", width = 10, height = 6, units = "in")

################

tmp <- t0_obj[["Lineage"]]@counts; tmp@x <- rep(1, length(tmp@x))
colsum_vec <- Matrix::colSums(tmp)
png("../../../../out/figures/Writeup4a/2022-02-09_lineage_histogram.png", 
    height = 900, width = 1500, res = 300, units = "px")
graphics::hist(colsum_vec, breaks = seq(min(colsum_vec)-.5, max(colsum_vec)+.5, by = 1),
               col = "gray", main = "Number of lineage barcodes for a cell",
               xlab = "Number of lineages", ylab = "Frequency")
graphics.off()

##################

t0_obj <- Seurat::SCTransform(t0_obj)
t0_obj <- Seurat::RunPCA(t0_obj, npcs = 50, verbose = FALSE)
set.seed(10)
t0_obj <- Seurat::RunUMAP(t0_obj, reduction = "pca", dims = 1:50)
t0_obj <- Seurat::FindNeighbors(t0_obj, reduction = "pca", dims = 1:50, verbose = FALSE)
t0_obj <- Seurat::FindClusters(t0_obj, resolution = 1, verbose = FALSE)

plot1 <- Seurat::DimPlot(t0_obj, 
                      reduction = "umap", 
                      group.by = "seurat_clusters", 
                      label = TRUE,
                      repel = TRUE)
plot2 <- Seurat::FeaturePlot(t0_obj, 
                          reduction = "umap", 
                          features = "nCount_RNA")
lineage_cutoff <- t0_obj$nCount_Lineage
lineage_cutoff[lineage_cutoff > 1] <- 2
t0_obj$Lineage_cutoff <- lineage_cutoff
plot3 <- Seurat::FeaturePlot(t0_obj, 
                             reduction = "umap", 
                             features = "Lineage_cutoff")
plot1 <- plot1 + plot2 + plot3
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_umap.png"),
                plot1, device = "png", width = 15, height = 4.5, units = "in")

###############

mito_idx <- grep("^MT-", rownames(t0_obj[["SCT"]]@scale.data))
rownames(t0_obj[["SCT"]]@scale.data)[mito_idx]
mito_l2 <- sqrt(Matrix::colSums(t0_obj[["SCT"]]@scale.data[mito_idx,]^2))
all_l2 <- sqrt(Matrix::colSums(t0_obj[["SCT"]]@scale.data^2))
ratio_vec <- mito_l2/all_l2
t0_obj$mt.ratio.sctransform <- ratio_vec
plot1 <- Seurat::FeaturePlot(t0_obj, 
                             reduction = "umap", 
                             features = c("percent.mt", "mt.ratio.sctransform"),
                             ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_umap_mitochondria.png"),
                plot1, device = "png", width = 10, height = 4.5, units = "in")

mito_idx <- grep("^MT-", rownames(t0_obj[["RNA"]]@counts))
mito_sum <- Matrix::colSums(t0_obj[["RNA"]]@counts[mito_idx,])
t0_obj$nCount_RNA_mito <- mito_sum

Seurat::Idents(t0_obj) <- rep("1", nrow(t0_obj[["RNA"]]))
plot1 <- Seurat::FeatureScatter(t0_obj, 
                                feature1 = "nCount_RNA_mito", feature2 = "nCount_RNA",
                                group.by = NULL)
plot1 <- plot1 + ggplot2::theme_gray() + Seurat::NoLegend() + ggplot2::coord_fixed()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_QC_scatterplot_mito.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

############

res <- test_barcode_assignment(seurat_obj = t0_obj)
t0_obj <- res$seurat_obj
lineage_preprocessing_outs <- res
lineage_preprocessing_outs$seurat_obj <- NULL

save(t0_obj, lineage_preprocessing_outs, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4a/2022-02-09_time0_preprocessed.RData")

plot1 <- plot_barcode_threshold(res$threshold_df)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_barcode_threshold.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")
plot1 <- plot_lineage_barcodecounts(t0_obj)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_lineage_count.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

