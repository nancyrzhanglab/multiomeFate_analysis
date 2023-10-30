rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6o/Writeup6o_all-data_lightweight.RData")

tmp <- all_data[["lsi"]]
tmp@assay.used <- "Saver"

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["lsi"]] <- NULL
all_data[["lsi"]] <- tmp

all_data[["ATAC"]] <- NULL
all_data[["RNA"]] <- NULL
all_data[["Lineage"]] <- NULL
all_data[["geneActivity"]] <- NULL
all_data[["umap"]] <- NULL
all_data[["atac.umap"]] <- NULL
all_data[["activity.umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["saverumap"]] <- NULL

