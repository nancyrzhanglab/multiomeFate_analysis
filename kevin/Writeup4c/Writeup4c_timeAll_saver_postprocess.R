rm(list=ls())
load("../../../../out/kevin/Writeup4c/Writeup4c_timeAll_sctransform_merging.RData")
load("../../../../out/kevin/Writeup4c/Writeup4c_timeAll_saver.RData")
library(Seurat)
library(SAVER)

quantile(saver_res$estimate, probs = seq(0.9, 1, length.out = 11))
tmp <- saver_res$estimate
tmp <- pmin(tmp, 10)
all_data[["Saver"]] <- Seurat::CreateAssayObject(counts = tmp)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data <- Seurat::ScaleData(all_data)
all_data[["Saver"]]@var.features <- rownames(tmp)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE,
                           reduction.name = "saverpca") 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 1:50,
                            reduction = "saverpca",
                            reduction.name = "saverumap")

plot1 <-Seurat::DimPlot(all_data, reduction = "saverumap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (SAVER),\n", length(all_data[["SCT"]]@var.features), " genes, using 50 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4c/Writeup4c_rna_umap_saver.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "saverumap",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E", 
                                               "CD44", "LOXL2", "ID3")),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4c/Writeup4c_rna_umap_saver_jackpot1.png"),
                plot1, device = "png", width = 12, height = 12, units = "in")

################################

rm(list=ls())
load("../../../../out/kevin/Writeup4c/Writeup4c_timeAll_saver.RData")
load("../../../../out/kevin/Writeup4c/Writeup4c_timeAll_peakmerging_simplified.RData")

library(Seurat)
library(Signac)
library(SAVER)

all_data[["Saver"]] <- Seurat::CreateAssayObject(counts = saver_res$estimate)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data <- Seurat::ScaleData(all_data)
all_data[["Saver"]]@var.features <- rownames(saver_res$estimate)
set.seed(10)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE,
                           reduction.name = "saverpca") 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 1:50,
                            reduction = "saverpca",
                            reduction.name = "saverumap")

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4c/Writeup4c_timeAll_ATAC-and-SAVER.RData")


