rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)

# see also https://github.com/nancyrzhanglab/multiomeFate_analysis/blob/emilia/emilia/task_identify_features_predictive_growth_and_lineage_specific/Run_chromVAR_Vierstra.R

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc-step2.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

## see https://stuartlab.org/signac/articles/motif_vignette.html
## https://stuartlab.org/signac/articles/data_structures.html#the-motif-class
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

main.chroms <- GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(GenomeInfoDb::seqnames(GenomicRanges::granges(all_data[["ATAC"]]))) %in% main.chroms)
all_data[["ATAC"]] <- subset(all_data[["ATAC"]], 
                             features = SeuratObject::Features(all_data[["ATAC"]])[keep.peaks])

all_data <- subset(all_data, dataset == "day10_CIS")

################################
# construct the counts
################################
peaks <- GenomicRanges::granges(all_data[["ATAC"]])
metadat <- all_data@meta.data
fm <- SeuratObject::LayerData(all_data, 
                              assay = "ATAC",
                              layer = "counts")

# Binarize 
fm@x <- rep(1, length(fm@x)) # why?

counts <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = fm), 
  rowRanges = peaks, 
  colData = metadat)

# Remove non-standard chromosomes
chr.exclude <- GenomeInfoDb::seqlevels(SummarizedExperiment::rowRanges(counts))[grep("random|chrM|Un|GL|KI", GenomeInfoDb::seqlevels(SummarizedExperiment::rowRanges(counts)))]
if(length(chr.exclude) > 0) {
  idy = grep(paste(chr.exclude, collapse="|"), SummarizedExperiment::rowRanges(counts))
  if(length(idy) > 0) {counts = counts[-idy, ]}
}

# Remove peaks with 0 fragments across all cells
remove_idx <- as.numeric(which(Matrix::rowSums(SummarizedExperiment::assay(counts)) == 0))
if(length(remove_idx) > 0) counts <- counts[-remove_idx, ]

################################
# run ChromVar
################################

# Calculate GC bias
counts <- chromVAR::addGCBias(
  counts, 
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
)

motif.positions <- motifmatchr::matchMotifs(
  pwms = pfm,
  subject = counts,
  out = "positions",
  genome = "hg38"
)

dev <- chromVAR::computeDeviations(
  object = counts, 
  annotations = motif.positions
)

mat.motifs <- chromVAR::deviationScores(dev)


