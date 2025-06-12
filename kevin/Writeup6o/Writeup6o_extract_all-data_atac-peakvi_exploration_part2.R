rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

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

print("Saving")
save(date_of_run, session_info,
     all_data,
     file = "../../../../out/kevin/Writeup6o/Writeup6o_extract_all-data_atac-peakvi_exploration.RData")


