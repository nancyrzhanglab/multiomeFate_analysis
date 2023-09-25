rm(list=ls())
library(Seurat); library(Signac)
library(GenomicRanges); library(GenomeInfoDb); library(IRanges)
library(JASPAR2020); library(TFBSTools); library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges); library(IRanges)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data) <- "ATAC"
all_data <- Signac::RunChromVAR(
  object = all_data,
  genome = "hg38"
)

Seurat::DefaultAssay(all_data) <- "Saver"
note <- paste0(note, " Then, adding chromVar scores were added to Writeup6m_all-data_add-motif.RData. ")

print("Saving")
save(date_of_run, session_info, 
     all_data, note,
     file = "../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")
