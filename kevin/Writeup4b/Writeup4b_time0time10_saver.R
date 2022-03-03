rm(list=ls())
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_exploration_alldata_onlyGEX.RData")
library(Seurat)
library(SAVER)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
saver_res <- SAVER::saver(x = mat, ncores = 1)
save(saver_res, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_saver.RData")

