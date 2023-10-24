rm(list=ls())
library(Seurat)
library(Signac)
library(SeuratDisk)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

all_data

print("Removing unnecessary modalities")
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["activityPCA"]] <- NULL
all_data[["activity.umap"]] <- NULL
all_data[["saverumap"]] <- NULL
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
all_data[["RNA"]] <- NULL
all_data[["Saver"]] <- NULL
all_data[["Lineage"]] <- NULL
all_data[["geneActivity"]] <- NULL

all_data
all_data[["ATAC"]]

print("Subsetting")
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_CIS", "week5_CIS"))] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)
all_data[["ATAC"]]@motifs <- NULL # see https://github.com/mojaveazure/seurat-disk/issues/15#issuecomment-1544286445
all_data[["ATAC"]]@positionEnrichment <- list()
all_data[["ATAC"]]@data <- Matrix::Matrix(0, nrow = 2, ncol = ncol(all_data), sparse = T)
all_data[["ATAC"]]@scale.data <- matrix(0, nrow = 2, ncol = ncol(all_data))

print("Saving")
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
SeuratDisk::SaveH5Seurat(all_data, filename = "../../../../out/kevin/Writeup6o/Writeup6o_all-data-atac_CIS.h5Seurat")
SeuratDisk::Convert("../../../../out/kevin/Writeup6o/Writeup6o_all-data-atac_CIS.h5Seurat", dest = "h5ad", misc = F)

print("Done! :)")