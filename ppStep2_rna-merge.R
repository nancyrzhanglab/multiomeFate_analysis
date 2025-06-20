rm(list=ls())
library(Seurat)
library(SeuratObject)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(GenomicRanges)
library(IRanges)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

create_seurat_object <- function(file_prefix, file_folder, file_suffix){
  tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = tmp[["Gene Expression"]],
    assay = "RNA"
  )
  
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj
}

############

# file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_prefix <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2024_08/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2024_08_21_GEXandLin_time0", "2024_08_21_GEXandLin_day10_CIS", 
                  "2024_08_21_GEXandLin_day10_COCL2", "2024_08_21_GEXandLin_day10_DABTRAM",
                  "2024_08_21_GEXandLin_week5_CIS", "2024_08_21_GEXandLin_week5_COCL2",
                  "2024_08_21_GEXandLin_week5_DABTRAM")
name_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM",
              "week5_CIS", "week5_COCL2", "week5_DABTRAM")
out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

seurat_list <- lapply(1:length(file_folders), function(i){
  file_folder <- file_folders[i]
  seurat_obj <- create_seurat_object(file_prefix = file_prefix,
                                     file_folder = file_folder,
                                     file_suffix = file_suffix)
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})
  
all_data_rna <- merge(seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]], 
                                              seurat_list[[4]], seurat_list[[5]], 
                                              seurat_list[[6]], seurat_list[[7]]), 
                      add.cell.ids = name_vec, 
                      project = "All_Data", merge.data = TRUE)

print("Finished processing RNA")
save(all_data_rna, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep2_rna-merge.RData"))
