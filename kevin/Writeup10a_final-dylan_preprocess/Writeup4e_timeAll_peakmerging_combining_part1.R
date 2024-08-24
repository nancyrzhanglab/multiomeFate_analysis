rm(list=ls())

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# also port in the SAVER
all_data2 <- all_data
load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_saver.RData")
all(colnames(all_data) == colnames(all_data2))
keep_vec <- rep(0, ncol(all_data2))
keep_vec[colnames(all_data2) %in% colnames(all_data)] <- 1
all_data2$keep <- keep_vec
all_data2 <- subset(all_data2, keep == 1)
all(colnames(all_data) == colnames(all_data2))

save(all_data2, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete_tmp.RData")

all_data2[["Saver"]] <- all_data[["Saver"]]
all_data2[["saverpca"]] <- all_data[["saverpca"]]
all_data2[["saverumap"]] <- all_data[["saverumap"]]

all_data <- all_data2
save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete_tmp.RData")


