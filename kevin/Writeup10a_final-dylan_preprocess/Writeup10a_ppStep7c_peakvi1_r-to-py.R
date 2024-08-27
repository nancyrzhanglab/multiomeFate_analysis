rm(list=ls())
library(Seurat)
library(Signac)
library(SeuratDisk)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc-step2.RData"))

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["Lineage"]] <- NULL
all_data[["RNA"]] <- NULL

treatment_vec <- c("All", "CIS", "COCL2", "DABTRAM")
all_data_full <- all_data

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  all_data <- all_data_full
  
  print("Subsetting")
  if(treatment != "all"){
    keep_vec <- rep(FALSE, length(Seurat::Cells(all_data)))
    keep_vec[which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))] <- TRUE
    all_data$keep <- keep_vec
    all_data <- subset(all_data, keep == TRUE)
  }
  all_data[["ATAC"]]@motifs <- NULL # see https://github.com/mojaveazure/seurat-disk/issues/15#issuecomment-1544286445
  all_data[["ATAC"]]@positionEnrichment <- list()
  all_data <- Seurat::DietSeurat(all_data, 
                                 layers = "counts")
  
  print("Saving")
  # https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
  SeuratDisk::SaveH5Seurat(all_data, 
                           filename = paste0(out_folder, "Writeup10a_ppStep7c_peakvi-prep_", treatment, ".h5Seurat"))
  SeuratDisk::Convert(paste0(out_folder, "Writeup10a_ppStep7c_peakvi-prep_", treatment, ".h5Seurat"), 
                      dest = "h5ad", 
                      misc = FALSE)
}

print("Done! :)")