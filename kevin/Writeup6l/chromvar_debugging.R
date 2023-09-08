rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(SummarizedExperiment)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract_lightweight.RData")

## from https://github.com/stuart-lab/signac/blob/2ad6c3c9c0c8dd31f7e1433b2efd5050d8606f27/R/motifs.R#L150
set.seed(10)
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["chromvar"]] <- NULL
all_data[["ATAC"]]@motifs <- NULL

# perhaps consider this one: https://compgenomr.github.io/book/motif-discovery.html
print("Apply TFBSTools::getMatrixSet")
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

print("Apply Signac::CreateMotifMatrix")
motif.matrix <- Signac::CreateMotifMatrix(
  features = GenomicRanges::granges(all_data[["ATAC"]]),
  pwm = pfm,
  genome = "hg38",
  use.counts = FALSE
)

print("Apply motifmatchr::matchMotifs")
motif.positions <- motifmatchr::matchMotifs(
  pwms = pfm,
  subject = GenomicRanges::granges(all_data[["ATAC"]]),
  out = "positions",
  genome = "hg38"
)

print("Apply Signac::CreateMotifObject")
motif <- Signac::CreateMotifObject(
  data = motif.matrix,
  positions = motif.positions,
  pwm = pfm
)

print("Apply Seurat::SetAssayData")
all_data[["ATAC"]] <- Seurat::SetAssayData(
  object = all_data[["ATAC"]],
  slot = 'motifs',
  new.data = motif
)

##################

motif.matrix <- GetMotifData(object = all_data, slot = "data")
peak.matrix <- GetAssayData(object = all_data, slot = "counts")

idx.keep <- rowSums(x = peak.matrix) > 0
peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
motif.matrix <- motif.matrix[idx.keep, , drop = FALSE]
peak.ranges <- granges(x = all_data)
peak.ranges <- peak.ranges[idx.keep]

chromvar.obj <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = peak.matrix),
  rowRanges = peak.ranges
)

chromvar.obj <- chromVAR::addGCBias(
  object = chromvar.obj,
  genome = "hg38"
)
head(chromvar.obj@rowRanges@elementMetadata)

row.data <- data.frame(SummarizedExperiment::rowData(x = chromvar.obj))
row.data[is.na(x = row.data)] <- 0
SummarizedExperiment::rowData(x = chromvar.obj) <- row.data

bg <- chromVAR::getBackgroundPeaks(
  object = chromvar.obj
)

######

dim(bg)
bg[1:5,1:5]

#############################

object <- chromvar.obj
bias <- SummarizedExperiment::rowData(object)$bias
niterations = 50
w = 0.1
bs = 50

fragments_per_peak <- chromVAR:::getFragmentsPerPeak(object)

intensity <- log10(fragments_per_peak)
norm_mat <- matrix(c(intensity, bias), ncol = 2, byrow = FALSE)

chol_cov_mat <- chol(cov(norm_mat))
trans_norm_mat <- t(forwardsolve(t(chol_cov_mat), t(norm_mat)))
apply(trans_norm_mat, 2, sd) # both 1's

bins1 <- seq(min(trans_norm_mat[, 1]), max(trans_norm_mat[, 1]), 
             length.out = bs)
bins2 <- seq(min(trans_norm_mat[, 2]), max(trans_norm_mat[, 2]), 
             length.out = bs)

bin_data <- do.call(rbind, lapply(seq_len(bs), 
                                  function(x) matrix(c(rep(bins1[x], bs), 
                                                       bins2), ncol = 2, 
                                                     byrow = FALSE)))
dim(bin_data)

bin_dist <- chromVAR:::euc_dist(bin_data)
bin_dist[1:5,1:5]
bin_p <- dnorm(bin_dist, 0, w)

bin_membership <- nabor::knn(bin_data, query = trans_norm_mat, k = 1)$nn.idx

bin_density <- chromVAR:::tabulate2(bin_membership, min_val = 1, max_val = bs^2)

background_peaks <- chromVAR:::bg_sample_helper(bin_membership - 1, bin_p, bin_density, 
                                                niterations)
dim(background_peaks)
background_peaks[1:5,1:5]
quantile(background_peaks)
dim(peak.matrix)
length(bin_membership)

x <- seq(1:10); x <- x/sum(x)
y <- chromVAR:::ProbSampleReplace(100000, x)

###################################

expectation <- chromVAR:::getFragmentsPerPeak(chromvar.obj) / chromVAR:::getTotalFragments(chromvar.obj)
zz <- chromVAR:::getTotalFragments(chromvar.obj)

