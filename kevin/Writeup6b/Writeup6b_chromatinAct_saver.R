rm(list=ls())
library(Seurat)
library(Signac)
library(SAVER)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

mat <- all_data[["geneActivity"]]@counts[Seurat::VariableFeatures(all_data, assay = "geneActivity"),]
print(dim(mat))
saver_res <- SAVER::saver(x = mat, ncores = 4)
save(saver_res, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_chromatinAct_saver.RData")
