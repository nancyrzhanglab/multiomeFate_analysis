library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(HDF5Array)

in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/LARRY_hematopoiesis/'
load(paste0(in_dir, "Writeup8_larry-dataset.RData"))
load(paste0(in_dir,'/Writeup8_larry-dataset_step2_lineage-plotting.RData'))
seurat_object[["RNA"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")
# seurat_object@reductions[["umap"]] <- NULL
# seurat_object@reductions[["SPRING"]] <- NULL
SaveH5Seurat(seurat_object, filename = paste0(in_dir, "Writeup8_larry-dataset_step2_lineage-plotting-spring.h5Seurat"))
Convert(paste0(in_dir, "Writeup8_larry-dataset_step2_lineage-plotting-spring.h5Seurat"), dest = "h5ad")

DimPlot(seurat_object, reduction = 'SPRING', split.by = 'Time.point')

seurat_object.meta <- seurat_object@meta.data
write.csv(seurat_object.meta, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/scripts/task5_cospar/kevin_larry_data_meta.csv')
