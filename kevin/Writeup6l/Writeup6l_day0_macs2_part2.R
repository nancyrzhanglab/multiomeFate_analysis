rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract.RData")
load("../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

fragments <- Signac::CreateFragmentObject("~/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time0/outs/atac_fragments.tsv.gz")

set.seed(10)
time_start <- Sys.time()
feature_mat <- Signac::FeatureMatrix(
  fragments = fragments,
  features = peaks,
  cells = colnames(all_data)
)
time_end <- Sys.time()

save(feature_mat, 
     date_of_run, session_info, time_start, time_end,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2_part2.RData")
