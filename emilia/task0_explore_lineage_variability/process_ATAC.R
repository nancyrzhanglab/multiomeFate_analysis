library(Seurat)
library(Signac)
library(IRanges)

args <- commandArgs(trailingOnly=TRUE)
treatment <- args[1]

# ==============================================================================
# Load data
# ==============================================================================

data <- readRDS('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/data/TBK1i_Multiome_StringentMT_InVivo_Run2_ATAC_Cancer_Harmony_eigs12.scMultiome_combined_object_Filtered.rds')
Seurat::DefaultAssay(data) <- "ATAC"

data$pca <- NULL
data$umap <- NULL

idx <- which(data$Sample == treatment)
data <- data[,idx]

# ==============================================================================
# Process data (ATAC only)
# ==============================================================================

data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunSVD(data)
data <- RunUMAP(data, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# DimPlot(data, reduction = "umap.atac", split.by = 'Treatment')

data_lsi <- data@reductions[["lsi"]]@cell.embeddings
write.csv(data_lsi, paste0('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/data/TBK1i_Multiome_StringentMT_InVivo_Run2_ATAC_Cancer_Harmony_eigs12.scMultiome_combined_object_Filtered_', treatment, '_lsi.csv'))

