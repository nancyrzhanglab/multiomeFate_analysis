rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)

load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data_subset) <- "geneActivity"
mat <- all_data_subset[["geneActivity"]]@counts[Seurat::VariableFeatures(all_data_subset),]
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
time_start <- Sys.time()
topic_res <- fastTopics::fit_topic_model(mat, k = K)
time_end <- Sys.time()

save(topic_res, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_fasttopics_DABTRAM-geneactivity.RData")


