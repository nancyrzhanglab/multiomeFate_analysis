rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract.RData")

Seurat::DefaultAssay(all_data) <- "ATAC"
all_data <- Signac::RunSVD(all_data)
rann_result <- RANN::nn2(all_data[["lsi"]]@cell.embeddings[,2:30])

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["ATAC"]] <- NULL

save(date_of_run, session_info, 
     all_data, rann_result,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_atac_NN.RData")