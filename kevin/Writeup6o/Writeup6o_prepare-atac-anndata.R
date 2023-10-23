rm(list=ls())
library(Seurat)
library(Signac)
library(SeuratDisk)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

all_data

print("Removing unnecessary modalities")
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["pca"]] <- NULL
all_data[["umap"]] <- NULL
all_data[["lsi"]] <- NULL
all_data[["atac.umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["chromvar"]] <- NULL
all_data[["RNA"]] <- NULL
all_data[["Saver"]] <- NULL

print("Subsetting")
all_data <- subset(all_data, dataset == "CIS")

print("Saving")
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
SeuratDisk::SaveH5Seurat(pbmc3k.final, filename = "../../../../out/kevin/Writeup6o/Writeup6o_all-data-atac_CIS.h5Seurat")
SeuratDisk::Convert("../../../../out/kevin/Writeup6o/Writeup6o_all-data-atac_CIS.h5Seurat", dest = "h5ad")

print("Done! :)")