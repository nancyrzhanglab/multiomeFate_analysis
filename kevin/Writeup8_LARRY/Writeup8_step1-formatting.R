rm(list=ls())
library(Seurat)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset.RData")

# remove all the cells not in the trajectory
seurat_object <- subset(seurat_object, trajectory == TRUE)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Subsetting out to only keep the cells in the trajectory.")

save(seurat_object,
     date_of_run, session_info, note,
     file = "~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_subset.RData")

print("Done! :)")
