rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6h/peak-gene-matching.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
metadata <- all_data@meta.data

treatment <- "DABTRAM"
keep_vec <- rep(F, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_DABTRAM", "week5_DABTRAM"))] <- T
all_data$keep <- keep_vec

all_data <- subset(all_data, keep == TRUE)
metadata_mat <- all_data@meta.data
print(all_data$dataset)

preprocess_res <- multiomeFate:::extract_relevant_peaks(peak_mapping_list = matching_list,
                                         seurat_obj = all_data,
                                         slot_atac = "ATAC",
                                         slot_rna = "Saver",
                                         verbose = 3)
preprocess_res <- multiomeFate:::preprocess_chromatin_peak(chr_peak_list = preprocess_res$chr_peak_list,
                                            rna_mat = preprocess_res$rna_mat,
                                            verbose = 3)

save(preprocess_res, tab_mat, metadata_mat,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM_supervised-pca_preprocess.RData")
