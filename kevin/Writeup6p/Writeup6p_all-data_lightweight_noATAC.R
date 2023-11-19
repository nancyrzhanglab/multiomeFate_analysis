rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data) <- "RNA"

all_data[["peakVI_CIS"]]@assay.used <- "RNA"
all_data[["peakVI_COCL2"]]@assay.used <- "RNA"
all_data[["peakVI_DABTRAM"]]@assay.used <- "RNA"
all_data[["ATAC"]] <- NULL

print("Saving")
save(date_of_run, session_info,
     all_data,
     file = "../../../../out/kevin/Writeup6p/Writeup6p_all-data_lightweight_noATAC.RData")

print("Done! :)")