rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
source("seurat_helpers.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_02_arc_time0", "2022_02_arc_time10_CIS", 
                  "2022_02_arc_time10_COCL2", "2022_02_arc_time10_DABTRAM")

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

time0 <- create_seurat_object(file_folders[1])
time10_cis <- create_seurat_object(file_folders[2])
time10_cocl2 <- create_seurat_object(file_folders[3])
time10_dabtram <- create_seurat_object(file_folders[4])

# https://satijalab.org/signac/articles/merging.html