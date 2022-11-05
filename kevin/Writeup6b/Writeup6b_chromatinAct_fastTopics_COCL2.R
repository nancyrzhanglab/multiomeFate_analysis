rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)
library(fastTopics)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))] <- 1
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == 1)

mat <- all_data[["geneActivity"]]@counts[Seurat::VariableFeatures(all_data, assay = "geneActivity"),]
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
time_start <- Sys.time()
topic_res <- fastTopics::fit_topic_model(mat, k = K)
time_end <- Sys.time()

save(topic_res, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_chromatinAct_fasttopics_COCL2.RData")


