rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_processed.RData")
tmp <- all_data

load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_exploration_alldata_onlyGEX.RData")
all_data[["ATAC"]] <- tmp[["ATAC"]]
covariates <- c("nucleosome_signal", "nucleosome_percentile",
                "nucleosome_group", "TSS.enrichment",
                "TSS.percentile", "high.tss", 
                "blacklist_fraction", "passed_filters",
                "peak_region_fragments", "pct_reads_in_peaks")
all_data@meta.data <- cbind(all_data@meta.data, tmp@meta.data[,covariates])

# process the ATAC data
Seurat::DefaultAssay(all_data) <- "ATAC"
set.seed(10)
all_data <- Signac::RunTFIDF(all_data)
all_data <- Signac::FindTopFeatures(all_data, min.cutoff = 'q0')
all_data <- Signac::RunSVD(all_data)

set.seed(10)
all_data <- Seurat::RunUMAP(object = all_data, 
                            reduction = 'lsi', 
                            assay = 'ATAC', 
                            reduction.name = 'atac.umap', 
                            reduction.key = 'atacUMAP_',
                            dims = 2:30)

##############################

plot1 <-Seurat::DimPlot(all_data, reduction = "umap",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\n", length(all_data[["RNA"]]@var.features), " genes, using 25 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/tmp.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

set.seed(10)
all_data <- Seurat::FindMultiModalNeighbors(
  all_data, reduction.list = list("pca", "lsi"), 
  dims.list = list(1:25, 2:30), modality.weight.name = "RNA.weight"
)
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            nn.name = "weighted.nn", 
                            reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_processed2.RData")

plot1 <-Seurat::DimPlot(all_data, reduction = "wnn.umap",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("WNN"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_wnn_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::VlnPlot(all_data, features = "RNA.weight", 
                         group.by = 'original_dataset', 
                         pt.size = 0.1) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_wnn_violin.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

