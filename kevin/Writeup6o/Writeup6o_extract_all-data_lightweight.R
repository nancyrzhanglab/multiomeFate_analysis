rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[sample(1:length(keep_vec), 100)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

print("Saving")
save(date_of_run, session_info,
     all_data,
     file = "../../../../out/kevin/Writeup6o/Writeup6o_all-data_lightweight.RData")
