rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/Writeup6b/Writeup6b_chromVar.RData")
Seurat::DefaultAssay(all_data) <- "Saver"
print("1")
all_data[["ATAC"]] <- NULL
print("2")
all_data[["geneActivity"]] <- NULL
print("3")
all_data[["pca"]] <- NULL
print("4")
all_data[["umap"]] <- NULL
print("5")
all_data[["saverpca"]] <- NULL
print("6")
all_data[["fasttopic_CIS"]] <- NULL
print("7")
all_data[["fasttopic_COCL2"]] <- NULL
print("8")
all_data[["fasttopic_DABTRAM"]] <- NULL
print("9")
all_data[["common_tcca"]] <- NULL
print("10")
all_data[["distinct1_tcca"]] <- NULL
print("11")
all_data[["distinct2_tcca"]] <- NULL
print("12")
all_data[["activityPCA"]] <- NULL

save(date_of_run, session_info, all_data,
     file = "../../out/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")