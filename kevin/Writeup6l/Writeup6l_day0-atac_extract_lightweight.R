rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract.RData")

keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[sample(1:length(keep_vec), 100)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

print("Saving")
save(date_of_run, session_info, 
     all_data, 
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract_lightweight.RData")