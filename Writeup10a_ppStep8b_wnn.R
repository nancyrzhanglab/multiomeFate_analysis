rm(list=ls())
library(Seurat)
library(Signac)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"

all_data <- data_loader(which_files = c("rna", "saver", "peakvi"))

set.seed(10)
all_data <- Seurat::FindMultiModalNeighbors(
  all_data, 
  reduction.list = list("Saver.pca", "peakVI.All"), 
  dims.list = list(1:30, 1:14), 
  modality.weight.name = "RNA.weight"
)

set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            nn.name = "weighted.nn", 
                            reduction.name = "wnn.umap", 
                            reduction.key = "wnnUMAP_")

all_data_wnn <- all_data[["wnn.umap"]]

date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(all_data_wnn, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_wnn.RData"))
