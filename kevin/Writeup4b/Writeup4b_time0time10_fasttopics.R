rm(list=ls())
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_exploration_alldata_onlyGEX.RData")
library(Seurat)
library(fastTopics)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
time_start <- Sys.time()
topic_res <- fastTopics::fit_topic_model(mat, k = K)
time_end <- Sys.time()

save(topic_res, all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_fasttopics_onlyGEX.RData")

