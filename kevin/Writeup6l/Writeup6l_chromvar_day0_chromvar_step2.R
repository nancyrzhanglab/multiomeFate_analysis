rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromvar_day0.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

print("Extracting motifs")
data.use <- Signac::GetMotifData(object = all_data,
                                 assay = "ATAC",
                                 slot = "pwm")
names(data.use) <- Signac::GetMotifData(
  object = all_data,
  assay = "ATAC",
  slot = "motif.names"
)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["ATAC"]] <- NULL

save(date_of_run, session_info, 
     all_data, data.use,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0_chromvar_lightweight_noATAC.RData")

