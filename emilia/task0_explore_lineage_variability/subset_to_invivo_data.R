library(Seurat)
library(Signac)

sc_obj <- readRDS('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/data/TBK1i_Multiome_StringentMT_InVivo_Run2_ATAC_Cancer_Harmony_eigs12.scMultiome_combined_object_Filtered.rds')

meta <- sc_obj@meta.data
write.csv(meta, '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task0_explore_lineage_variability/data/TBK1i_Multiome_StringentMT_InVivo_Run2_ATAC_Cancer_Harmony_eigs12.scMultiome_combined_object_Filtered_meta.csv')

# data <- as(sc_obj@assays[["ChromVAR_Signatures"]]@data, "sparseMatrix") 
# saveRDS(data,'TBK1i_Multiome_StringentMT_InVivo_Run2_ATAC_Cancer_Harmony_eigs12.scMultiome_combined_object_Filtered_ChromVAR_Signatures.rds')