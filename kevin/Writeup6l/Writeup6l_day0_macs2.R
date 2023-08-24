rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# fix the paths
all_data[["ATAC"]]@fragments[[1]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time0/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[2]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_CIS/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[3]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_COCL2/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[4]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_DABTRAM/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[5]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_CIS/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[6]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_COCL2/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[7]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_DABTRAM/outs/atac_fragments.tsv.gz"

# https://stuartlab.org/signac/articles/peak_calling.html
# https://stuartlab.org/signac/reference/callpeaks

# be sure to have activated the python virtual environment: source venvMacs2/bin/activate
set.seed(10)
time_start <- Sys.time()
peaks <- Signac::CallPeaks(
  object = all_data[["ATAC"]]
)
time_end <- Sys.time()

save(peaks, 
     date_of_run, session_info, time_start, time_end,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2.RData")

########

length(peaks)
head(peaks)
levels(peaks@seqnames)
table(peaks@seqnames)
