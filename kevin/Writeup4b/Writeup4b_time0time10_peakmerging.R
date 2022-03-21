rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
source("seurat_helpers.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_02_arc_time0", "2022_02_arc_time10_CIS", 
                  "2022_02_arc_time10_COCL2", "2022_02_arc_time10_DABTRAM")

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

time0 <- create_seurat_object(file_folders[1])
time10_cis <- create_seurat_object(file_folders[2])
time10_cocl2 <- create_seurat_object(file_folders[3])
time10_dabtram <- create_seurat_object(file_folders[4])

time0b <- time0; time0b[["ATAC"]] <- NULL
time10b_cis <- time10_cis; time10b_cis[["ATAC"]] <- NULL
time10b_cocl2 <- time10_cocl2; time10b_cocl2[["ATAC"]] <- NULL
time10b_dabtram <- time10_dabtram; time10b_dabtram[["ATAC"]] <- NULL

time0b$dataset <- "time0"
time10b_cis$dataset <- "time10_cis"
time10b_cocl2$dataset <- "time10_cocl2"
time10b_dabtram$dataset <- "time10_dabtram"

all_data <- merge(time0b, y = c(time10b_cis, time10b_cocl2, time10b_dabtram), 
                  add.cell.ids = c("time0", "time10_cis", "time10_cocl2", "time10_dabtram"), 
                  project = "All_Data", merge.data = T)

#######################

# https://satijalab.org/signac/articles/merging.html

# read in peak sets
peaks_time0 <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time0/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks_time10_cis <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_CIS/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks_time10_cocl2 <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_COCL2/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks_time10_dabtram <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_DABTRAM/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr_time0 <- GenomicRanges::makeGRangesFromDataFrame(peaks_time0)
gr_time10_cis <- GenomicRanges::makeGRangesFromDataFrame(peaks_time10_cis)
gr_time10_cocl2 <- GenomicRanges::makeGRangesFromDataFrame(peaks_time10_cocl2)
gr_time10_dabtram <- GenomicRanges::makeGRangesFromDataFrame(peaks_time10_dabtram)

# Create a unified set of peaks to quantify in each dataset
combined_peaks <- reduce(x = c(gr_time0, gr_time10_cis, gr_time10_cocl2, gr_time10_dabtram))

# Filter out bad peaks based on length
peak_widths <- IRanges::width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths < 10000 & peak_widths > 20]

# create fragment objects
frags_time0 <- Signac::CreateFragmentObject(
  path = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time0/outs/atac_fragments.tsv.gz",
  cells = colnames(time0)
)
frags_time10_cis <- Signac::CreateFragmentObject(
  path = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_CIS/outs/atac_fragments.tsv.gz",
  cells = colnames(time10_cis)
)
frags_time10_cocl2 <- Signac::CreateFragmentObject(
  path = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_COCL2/outs/atac_fragments.tsv.gz",
  cells = colnames(time10_cocl2)
)
frags_time10_dabtram <- Signac::CreateFragmentObject(
  path = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_DABTRAM/outs/atac_fragments.tsv.gz",
  cells = colnames(time10_dabtram)
)

time0_atac_counts <- Signac::FeatureMatrix(
  fragments = frags_time0,
  features = combined_peaks,
  cells = colnames(time0)
)
time10_cis_atac_counts <- Signac::FeatureMatrix(
  fragments = frags_time10_cis,
  features = combined_peaks,
  cells = colnames(time10_cis)
)
time10_cocl2_atac_counts <- Signac::FeatureMatrix(
  fragments = frags_time10_cocl2,
  features = combined_peaks,
  cells = colnames(time10_cocl2)
)
time10_dabtram_atac_counts <- Signac::FeatureMatrix(
  fragments = frags_time10_dabtram,
  features = combined_peaks,
  cells = colnames(time10_dabtram)
)

time0_atac_assay <- Signac::CreateChromatinAssay(counts = time0_atac_counts,
                                                 fragments = frags_time0,
                                                 sep = c("-", "-"),
                                                 annotation = annotation)
time0_atac <- Seurat::CreateSeuratObject(time0_atac_assay, assay = "ATAC")

time10_cis_atac_assay <- Signac::CreateChromatinAssay(counts = time10_cis_atac_counts,
                                                      fragments = frags_time10_cis,
                                                      sep = c("-", "-"),
                                                      annotation = annotation)
time10_cis_atac <- Seurat::CreateSeuratObject(time10_cis_atac_assay, assay = "ATAC")

time10_cocl2_atac_assay <- Signac::CreateChromatinAssay(counts = time10_cocl2_atac_counts,
                                                        fragments = frags_time10_cocl2,
                                                        sep = c("-", "-"),
                                                        annotation = annotation)
time10_cocl2_atac <- Seurat::CreateSeuratObject(time10_cocl2_atac_assay, assay = "ATAC")

time10_dabtram_atac_assay <- Signac::CreateChromatinAssay(counts = time10_dabtram_atac_counts,
                                                          fragments = frags_time10_dabtram,
                                                          sep = c("-", "-"),
                                                          annotation = annotation)
time10_dabtram_atac <- Seurat::CreateSeuratObject(time10_dabtram_atac_assay, assay = "ATAC")

time0_atac$dataset <- "time0"
time10_cis_atac$dataset <- "time10_cis"
time10_cocl2_atac$dataset <- "time10_cocl2"
time10_dabtram_atac$dataset <- "time10_dabtram"

all_data_atac <- merge(time0_atac, y = c(time10_cis_atac, time10_cocl2_atac, time10_dabtram_atac), 
                       add.cell.ids = c("time0", "time10_cis", "time10_cocl2", "time10_dabtram"), 
                       project = "All_Data", merge.data = T)

save(all_data, all_data_atac, 
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_peakmerging.RData")
