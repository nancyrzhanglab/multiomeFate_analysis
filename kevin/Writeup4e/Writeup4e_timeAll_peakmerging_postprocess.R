rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")


## RNA quality
Seurat::DefaultAssay(all_data) <- "RNA"
plot1 <- Seurat::VlnPlot(all_data, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
                         group.by = "dataset",
                         pt.size = 0.01,
                         ncol = 4)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_rna_QC.png"),
                plot1, device = "png", width = 20, height = 6, units = "in")

## ATAC quality
Seurat::DefaultAssay(all_data) <- "ATAC"
plot1 <- Signac::TSSPlot(all_data, group.by = 'high.tss') + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_tssEnrichment.png"),
                plot1, device = "png", width = 8, height = 4, units = "in")

plot1 <- Seurat::VlnPlot(
  object = all_data,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
  group.by = "dataset", 
  pt.size = 0.1,
  ncol = 5
)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_atac_QC.png"),
                plot1, device = "png", width = 20, height = 6, units = "in")

###### RNA UMAP

Seurat::DefaultAssay(all_data) <- "RNA"
plot1 <-Seurat::DimPlot(all_data, reduction = "umap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\n", length(all_data[["RNA"]]@var.features), " genes, using 50 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_rna_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

Seurat::DefaultAssay(all_data) <- "RNA"
plot1 <-Seurat::DimPlot(all_data, reduction = "umap",
                        group.by = "Phase", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\nCell-cycle phase"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_rna_umap-cellcycle.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

Seurat::DefaultAssay(all_data) <- "RNA"
plot1 <-Seurat::FeaturePlot(all_data, features = c("S.Score", "G2M.Score"))
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\nCell-cycle score"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_rna_umap-cellcycle-score.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")


plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E", 
                                               "CD44", "LOXL2", "ID3")),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_rna_umap_jackpot1.png"),
                plot1, device = "png", width = 12, height = 12, units = "in")

###### ATAC UMAP

Seurat::DefaultAssay(all_data) <- "ATAC"
plot1 <- Seurat::DimPlot(all_data, 
                         reduction = "atac.umap", 
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (TF-IDF),\nusing 49 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_atac_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###############################

Seurat::DefaultAssay(all_data) <- "geneActivity"
plot1 <- Seurat::DimPlot(all_data, 
                         reduction = "activity.umap", 
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (Gene activity),\nusing 49 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_atac-activity_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


