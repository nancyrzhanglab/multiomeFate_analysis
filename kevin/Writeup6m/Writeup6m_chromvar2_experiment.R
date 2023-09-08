rm(list=ls())
library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract_lightweight.RData")

set.seed(10)
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["chromvar"]] <- NULL
all_data[["ATAC"]]@motifs <- NULL

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

row.data <- data.frame(SummarizedExperiment::rowData(x = chromvar.obj))
row.data[is.na(x = row.data)] <- 0
SummarizedExperiment::rowData(x = chromvar.obj) <- row.data

bg <- chromVAR::getBackgroundPeaks(
  object = chromvar.obj,
  niterations = 200
)

##########################

tf <- "MA1132.1"
tf_idx <- which(colnames(motif.matrix) == tf)
peak_set <- which(motif.matrix[,tf_idx])

counts_mat <- chromvar.obj@assays@data$counts
tf_count <- length(peak_set)

tf_vec <- Matrix::sparseMatrix(j = peak_set,
                               i = rep(1, tf_count),
                               x = 1,
                               dims = c(1,
                                        nrow(counts_mat)))
observed <- as.vector(tf_vec %*% counts_mat)

background_peaks <- bg
niterations <- ncol(background_peaks)
sample_mat <- Matrix::sparseMatrix(j = as.vector(background_peaks[peak_set,
                                                                  seq_len(niterations)]),
                                   i = rep(seq_len(niterations), each = tf_count),
                                   x = 1,
                                   dims = c(niterations, nrow(counts_mat)))
sampled <- as.matrix(sample_mat %*% counts_mat)
sampled_mean <- Matrix::colMeans(sampled)
sampled_sd <- apply(sampled, 2, stats::sd)

z_score <- sapply(1:length(observed), function(i){
  (observed[i] - sampled_mean[i])/sampled_sd[i]
})

