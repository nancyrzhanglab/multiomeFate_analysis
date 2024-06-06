library(rhdf5)
library(SeuratDisk)

sceasy_seurat <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")
sceasy_seurat
file_name <- "/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/DABTRAM_data_MultiTimeClone_Later_FullSpace0_t*day0*day10*week5_adata_with_transition_map.h5ad"

Convert(file_name, dest = "h5seurat", overwrite = TRUE)
pbmc3k <- LoadH5Seurat("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/DABTRAM_data_MultiTimeClone_Later_FullSpace0_t*day0*day10*week5_adata_with_transition_map.h5seurat")
