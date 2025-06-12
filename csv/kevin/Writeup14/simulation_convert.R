library(Seurat)
library(SeuratDisk)
setwd("/home/mnt/weka/nzh/nzhanglab/project/Multiome_fate/out/kevin/Writeup14/")
load("Writeup14_plastic-setting_seurat_CoSpar_prepared.RData")
all_data[["Saver"]] <- as(object = all_data[["Saver"]], Class = "Assay")

all_data<- as.SingleCellExperiment(all_data)
library(zellkonverter)
out_path <- tempfile(pattern = ".h5ad")
writeH5AD(all_data, file = "Writeup14_plastic-setting_seurat_CoSpar_prepared.h5ad")
all_data@assays$RNA$scale.data ->scale.data
write.csv(scale.data, file ="plastic-setting_scale_data.csv",row.names=T,quote=F)

feature_mat = all_data[["fasttopic.COCL2"]]@cell.embeddings
write.csv(feature_mat, file ="plastic-setting_embeddings.csv",row.names=T,quote=F)



table_mat = table(all_data@meta.data$assigned_lineage,all_data@meta.data$dataset)
quantile(table_mat[,"week5_COCL2"])
table_mat = as.data.frame(table_mat)
table_mat$factor_lineage = as.character("High")
table_mat[table_mat$Freq<=92,]$factor_lineage <-  c("Low")
table_mat[table_mat$Var2==c("day10_COCL2"),]$factor_lineage <-  c("day10")
table_mat=subset(table_mat,table_mat$Var2==c("week5_COCL2"))
metadata<-all_data@meta.data
colnames(table_mat)=c("assigned_lineage","dataset","Freq","factor_lineage")
metadata = merge(metadata,table_mat[,c("assigned_lineage","factor_lineage")],by="assigned_lineage")
metadata[metadata$dataset==c("day10_COCL2"),]$factor_lineage <-  c("day10")
row.names(metadata)=row.names(all_data@meta.data)
write.csv(metadata,"plastic-setting_meta.csv",row.names=T,quote=F)

load("Writeup14_priming-setting_seurat_CoSpar_prepared.RData")
all_data@assays$RNA$scale.data ->scale.data
write.csv(scale.data, file ="priming-setting_scale_data.csv",row.names=T,quote=F)

feature_mat = all_data[["fasttopic.COCL2"]]@cell.embeddings
write.csv(feature_mat, file ="priming-setting_embeddings.csv",row.names=T,quote=F)

table_mat = table(all_data@meta.data$assigned_lineage,all_data@meta.data$dataset)
quantile(table_mat[,"week5_COCL2"])
table_mat = as.data.frame(table_mat)
table_mat$factor_lineage = as.character("High")
table_mat[table_mat$Freq<=107,]$factor_lineage <-  c("Low")
table_mat[table_mat$Var2==c("day10_COCL2"),]$factor_lineage <-  c("day10")
table_mat=subset(table_mat,table_mat$Var2==c("week5_COCL2"))
metadata<-all_data@meta.data
colnames(table_mat)=c("assigned_lineage","dataset","Freq","factor_lineage")
metadata = merge(metadata,table_mat[,c("assigned_lineage","factor_lineage")],by="assigned_lineage")
metadata[metadata$dataset==c("day10_COCL2"),]$factor_lineage <-  c("day10")
row.names(metadata)=row.names(all_data@meta.data)

write.csv(metadata,"priming-setting_meta.csv",row.names=T,quote=F)

all_data[["Saver"]] <- as(object = all_data[["Saver"]], Class = "Assay")

all_data<- as.SingleCellExperiment(all_data)
library(zellkonverter)
out_path <- tempfile(pattern = ".h5ad")
writeH5AD(all_data, file = "Writeup14_priming-setting_seurat_CoSpar_prepared.h5ad")


### prepare for PCA data ###
all_data.seurat <- as.Seurat(all_data)
DefaultAssay(object = all_data.seurat) <- "Saver"
all_data.seurat <- NormalizeData(all_data.seurat)
all_data.seurat <- FindVariableFeatures(object = all_data.seurat)
all_data.seurat <- ScaleData(object = all_data.seurat)
all_data.seurat <- RunPCA(object = all_data.seurat,npcs = 30)
mat_pca<- as.data.frame(as.matrix(all_data.seurat@reductions[["pca"]]@cell.embeddings))
saver_pca<- as.data.frame(as.matrix(all_data.seurat@reductions[["SAVER.PCA"]]@cell.embeddings))

write.csv(mat_pca, file = "pca30_priming-setting.csv",row.names=T,quote=F)
write.csv(saver_pca, file = "pca50_priming-setting.csv",row.names=T,quote=F)

load("Writeup14_plastic-setting_seurat_CoSpar_prepared.RData")
DefaultAssay(object = all_data) <- "Saver"
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(object = all_data)
all_data <- ScaleData(object = all_data)
all_data <- RunPCA(object = all_data,npcs = 30)
mat_pca<- as.data.frame(as.matrix(all_data@reductions[["pca"]]@cell.embeddings))
saver_pca<- as.data.frame(as.matrix(all_data@reductions[["Saver.pca"]]@cell.embeddings))

write.csv(mat_pca, file = "pca30_plastic-setting.csv",row.names=T,quote=F)
write.csv(saver_pca, file = "pca50_plastic-setting.csv",row.names=T,quote=F)

