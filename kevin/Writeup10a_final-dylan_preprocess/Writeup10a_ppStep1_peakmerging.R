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

# file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_prefix <- "/scratch/nzh/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_05_19_arc_time0", "2022_05_19_arc_time10_CIS", 
                  "2022_05_19_arc_time10_COCL2", "2022_05_19_arc_time10_DABTRAM",
                  "2022_05_19_arc_week5_CIS", "2022_05_19_arc_week5_COCL2",
                  "2022_05_19_arc_week5_DABTRAM")
name_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM",
              "week5_CIS", "week5_COCL2", "week5_DABTRAM")
out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

# https://satijalab.org/signac/articles/merging.html
# read in peak sets
print("Reading peak set")
peaks_list <- lapply(file_folders, function(file_folder){
  read.table(
    file = paste0(file_prefix, file_folder, "/outs/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
})

print("Creating genomic ranges")
# convert to genomic ranges
gr_list <- lapply(peaks_list, function(peaks_data){
  GenomicRanges::makeGRangesFromDataFrame(peaks_data)
})

print("Creating unified peak set")
# Create a unified set of peaks to quantify in each dataset
combined_peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], 
                               gr_list[[3]], gr_list[[4]], 
                               gr_list[[5]], gr_list[[6]], 
                               gr_list[[7]]))

# Filter out bad peaks based on length
peak_widths <- IRanges::width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths < 10000 & peak_widths > 20]

print("Creating frag objects")
# create fragment objects
frags_list <- lapply(1:length(file_folders), function(i){
  file_folder <- file_folders[i]
  Signac::CreateFragmentObject(
    path = paste0(file_prefix, file_folder, "/outs/atac_fragments.tsv.gz"),
    cells = colnames(seurat_list[[i]])
  )
})

print("Constructing count matrices")
atac_count_list <- lapply(1:length(file_folders), function(i){
  Signac::FeatureMatrix(
    fragments = frags_list[[i]],
    features = combined_peaks,
    cells = colnames(seurat_list[[i]])
  )
})

print("Constructing Seurat objects")
seurat_atac_list <- lapply(1:length(file_folders), function(i){
  assay_obj <- Signac::CreateChromatinAssay(counts = atac_count_list[[i]],
                                            fragments = frags_list[[i]],
                                            sep = c("-", "-"),
                                            annotation = annotation)
  seurat_obj <- Seurat::CreateSeuratObject(assay_obj, assay = "ATAC")
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})

print("Merging ATAC peaks")
all_data_atac <- merge(seurat_atac_list[[1]], y = c(seurat_atac_list[[2]], seurat_atac_list[[3]], 
                                                    seurat_atac_list[[4]], seurat_atac_list[[5]], 
                                                    seurat_atac_list[[6]], seurat_atac_list[[7]]), 
                       add.cell.ids = name_vec, 
                       project = "All_Data", merge.data = T)

print("Finished processing ATAC")
save(all_data_atac, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep1_peakmerging.RData"))
