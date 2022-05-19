rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folder <- "2022_04_28_GEXandLin_time0"
tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))

names(tmp)