rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)
load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))] <- 1
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == 1)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
time_start <- Sys.time()
topic_res <- fastTopics::fit_topic_model(mat, k = K)
time_end <- Sys.time()

save(topic_res, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_fasttopics_COCL2.RData")


