rm(list=ls())
load("../../../../out/kevin/Writeup4c/Writeup4c_timeAll_sctransform_merging.RData")
library(Seurat)
library(SAVER)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "SCT"),]
print(dim(mat))
saver_res <- SAVER::saver(x = mat, ncores = 4)
save(saver_res, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4c/Writeup4c_timeAll_saver.RData")

