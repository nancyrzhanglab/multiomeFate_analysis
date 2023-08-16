rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")
Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["ATAC"]] <- NULL
all_data[["geneActivity"]] <- NULL
all_data[["pca"]] <- NULL
all_data[["umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL

save(date_of_run, session_info, all_data,
     file = "../../out/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")