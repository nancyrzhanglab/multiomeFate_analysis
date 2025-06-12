# this simulation is for similar means but the different ranges

rm(list=ls())
library(Seurat)
library(multiomeFate)

load("~/project/Multiome_fate/out/kevin/Writeup9b/Writeup9b_simulation_day10-COCL2.RData")

set.seed(10)
rna_mat <- SeuratObject::LayerData(all_data, 
                                   layer = "counts", 
                                   assay = "RNA",
                                   features = Seurat::VariableFeatures(all_data))
saver_res <- SAVER::saver(x = rna_mat, ncores = 4)

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(saver_res,
     date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup9b/Writeup9b_day10-COCL2_saver.RData")

print("Done! :)")

