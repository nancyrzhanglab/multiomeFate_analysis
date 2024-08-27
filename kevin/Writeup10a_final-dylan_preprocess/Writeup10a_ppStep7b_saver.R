rm(list=ls())
library(Seurat)
library(Signac)
library(SAVER)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc-step2.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- SeuratObject::LayerData(all_data,
                               assay = "RNA",
                               layer = "counts",
                               features = Seurat::VariableFeatures(all_data, assay = "RNA"))
saver_res <- SAVER::saver(x = mat, ncores = 4)

print("Saving")
save(saver_res, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep7b_saver.RData"))

print("Done! :)")
