rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
mat <- Matrix::t(mat)

K <- 50
set.seed(10)
time_start <- Sys.time()
topic_res <- fastTopics::fit_topic_model(mat, k = K)
time_end <- Sys.time()

save(topic_res, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_RNA_fasttopics_all.RData")


