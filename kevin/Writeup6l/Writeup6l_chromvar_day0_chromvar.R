rm(list=ls())
library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)
Seurat::DefaultAssay(all_data) <- "ATAC"

## see https://stuartlab.org/signac/articles/motif_vignette.html
## https://stuartlab.org/signac/articles/data_structures.html#the-motif-class
print("Apply TFBSTools::getMatrixSet")
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# https://github.com/stuart-lab/signac/issues/486
print("Apply subset")
main.chroms <- GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(GenomeInfoDb::seqnames(GenomicRanges::granges(all_data[["ATAC"]]))) %in% main.chroms)
all_data[["ATAC"]] <- subset(all_data[["ATAC"]], features = rownames(all_data[["ATAC"]])[keep.peaks])

## see https://github.com/stuart-lab/signac/blob/master/R/motifs.R
print("Apply Signac::CreateMotifMatrix")
motif.matrix <- Signac::CreateMotifMatrix(
  features = GenomicRanges::granges(all_data[["ATAC"]]),
  pwm = pfm,
  genome = "hg38",
  use.counts = FALSE
)

save(date_of_run, session_info, 
     all_data, motif.matrix,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_chromvar_day0.RData")

###################################

print("Apply motifmatchr::matchMotifs")
motif.positions <- motifmatchr::matchMotifs(
  pwms = pfm,
  subject = GenomicRanges::granges(all_data[["ATAC"]]),
  out = "positions",
  genome = "hg38"
)

save(date_of_run, session_info, 
     all_data, motif.matrix, motif.positions,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_chromvar_day0.RData")

###################################

print("Apply Signac::CreateMotifObject")
motif <- Signac::CreateMotifObject(
  data = motif.matrix,
  positions = motif.positions,
  pwm = pfm
)

save(date_of_run, session_info, 
     all_data, motif.matrix, motif.positions, motif,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_chromvar_day0.RData")

###################################

print("Apply Seurat::SetAssayData")
all_data[["ATAC"]] <- Seurat::SetAssayData(
  object = all_data[["ATAC"]],
  slot = 'motifs',
  new.data = motif
)

print("Apply Signac::RunChromVAR")
all_data <- Signac::RunChromVAR(
  object = all_data,
  genome = "hg38"
)

save(date_of_run, session_info, 
     all_data, motif.matrix, motif.positions, motif,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_chromvar_day0.RData")

###################################

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["ATAC"]] <- NULL

save(date_of_run, session_info, 
     all_data, 
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0_chromvar_lightweight_noATAC.RData")

