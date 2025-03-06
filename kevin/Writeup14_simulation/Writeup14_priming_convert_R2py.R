# from https://github.com/nancyrzhanglab/multiomeFate_analysis/blob/kevin/csv/kevin/Writeup14/simulation_convert.R
rm(list=ls())

library(Seurat)
library(SeuratDisk)
library(zellkonverter)

# start of Sijia's code
setwd("/home/stat/nzh/team/kevinl1/nzhanglab/project/Multiome_fate/out/kevin/Writeup14/")

# write scale data
load("Writeup14_priming-setting_v2_seurat_CoSpar_prepared.RData")
scale.data <- all_data@assays$RNA$scale.data 
write.csv(scale.data, file ="priming-setting_v2_scale_data.csv", row.names = TRUE, quote = FALSE)

# write embedding
feature_mat <- all_data[["fasttopic.COCL2"]]@cell.embeddings
write.csv(feature_mat, file ="priming-setting_v2_embeddings.csv",row.names = TRUE,quote = FALSE)

# write tables
table_mat <- table(all_data@meta.data$assigned_lineage,all_data@meta.data$dataset)
cutoff <- sort(table_mat[,"week5_COCL2"], decreasing = TRUE)[10]
table_mat <- as.data.frame(table_mat)
table_mat$factor_lineage <- as.character("High")
table_mat[table_mat$Freq >= cutoff,]$factor_lineage <- c("Low")

table_mat[table_mat$Var2 == c("day10_COCL2"),]$factor_lineage <- c("day10")
table_mat <- subset(table_mat, table_mat$Var2 == c("week5_COCL2"))
metadata <- all_data@meta.data
colnames(table_mat) <- c("assigned_lineage","dataset","Freq","factor_lineage")
metadata <- merge(metadata,table_mat[,c("assigned_lineage","factor_lineage")],by="assigned_lineage")
metadata[metadata$dataset==c("day10_COCL2"),]$factor_lineage <- c("day10")
row.names(metadata) <- row.names(all_data@meta.data)
write.csv(metadata, "priming-setting_v2_meta.csv",row.names=T,quote=F)

# convert entire dataset into h5ad
all_data[["Saver"]] <- as(object = all_data[["Saver"]], Class = "Assay")
all_data_sce <- as.SingleCellExperiment(all_data)
out_path <- tempfile(pattern = ".h5ad")
zellkonverter::writeH5AD(all_data_sce, file = "Writeup14_priming-setting_v2_seurat_CoSpar_prepared.h5ad")