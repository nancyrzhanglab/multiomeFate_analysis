rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")
source("barcode_estimation.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

library_vec <- all_data$nCount_RNA
lin_mat <- all_data[["Lineage"]]@counts
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0, ]
dataset_vec <- factor(all_data$dataset)

keep_idx <- which(Matrix::colSums(lin_mat) > 0)
library_vec <- library_vec[keep_idx]
lin_mat <- lin_mat[,keep_idx]
dataset_vec <- dataset_vec[keep_idx]

barcode_res <- barcode_estimation(library_vec = library_vec,
                                  lin_mat = lin_mat,
                                  dataset_vec = dataset_vec,
                                  bool_shortcut = T,
                                  max_iter = 10,
                                  verbose = 2)

save(barcode_res, library_vec, lin_mat, dataset_vec,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4f/Writeup4f_barcode_estimation.RData")