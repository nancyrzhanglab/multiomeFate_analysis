rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("gene_list.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

## see https://stuartlab.org/signac/articles/motif_vignette.html
## https://stuartlab.org/signac/articles/data_structures.html#the-motif-class
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# https://github.com/stuart-lab/signac/issues/486
main.chroms <- GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(GenomeInfoDb::seqnames(GenomicRanges::granges(all_data[["ATAC"]]))) %in% main.chroms)
all_data[["ATAC"]] <- subset(all_data[["ATAC"]], features = rownames(all_data[["ATAC"]])[keep.peaks])

# annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86) # takes a while
# GenomeInfoDb::seqlevelsStyle(annotations) <- "NCBI" # chromosome 1 is called "chr1"
# ## see https://github.com/Bioconductor/GenomeInfoDb/issues/14 ?
# annotations
# chrom_assay <- all_data[["ATAC"]]
# Signac::Annotation(chrom_assay) <- annotations
# all_data[["ATAC"]] <- chrom_assay

Seurat::DefaultAssay(all_data) <- "ATAC"
# all_data <- Signac::AddMotifs(
#   object = all_data,
#   genome = "hg38",
#   pfm = pfm
# )

## see https://github.com/stuart-lab/signac/blob/master/R/motifs.R
motif.matrix <- Signac::CreateMotifMatrix(
  features = GenomicRanges::granges(all_data[["ATAC"]]),
  pwm = pfm,
  genome = "hg38",
  use.counts = FALSE
)

motif.positions <- motifmatchr::matchMotifs(
  pwms = pfm,
  subject = GenomicRanges::granges(all_data[["ATAC"]]),
  out = "positions",
  genome = "hg38"
)

motif <- Signac::CreateMotifObject(
  data = motif.matrix,
  positions = motif.positions,
  pwm = pfm
)

all_data[["ATAC"]] <- Seurat::SetAssayData(
  object = all_data[["ATAC"]],
  slot = 'motifs',
  new.data = motif
)

all_data <- Signac::RunChromVAR(
  object = all_data,
  genome = "hg38"
)