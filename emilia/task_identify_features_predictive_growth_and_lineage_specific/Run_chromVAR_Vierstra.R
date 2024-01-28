library(Seurat)
library(GenomicRanges)
library(chromVAR)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)

out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
# ==============================================================================
# Read data
# ==============================================================================
sc <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/day10_COCL2_processed.rds")
pwms <- readRDS("/Users/emiliac/Dropbox/Thesis/resources/motifs/Vierstra_Archetype_Motifs_v2.1.rds")

Seurat::DefaultAssay(sc) <- "ATAC"

peaks <- data.frame(
  seqnames = sapply(strsplit(rownames(sc), split="\\:|-"), function(x) x[[1]]),
  start = sapply(strsplit(rownames(sc), split="\\:|-"), function(x) x[[2]]),
  end = sapply(strsplit(rownames(sc), split="\\:|-"), function(x) x[[3]]))
peaks <- makeGRangesFromDataFrame(peaks)

metadat <- sc@meta.data
print(all(colnames(sc) == rownames(metadat)))

fm <- sc[['ATAC']]@counts

# ==============================================================================
# Cleaning data
# ==============================================================================
# Binarize #
fm <- ((fm > 0) + 0) # why?

# Convert to summarized experiment
counts <- SummarizedExperiment(
  assays = list(counts = fm), 
  rowRanges = peaks, 
  colData = metadat)

# Remove non-standard chromosomes
chr.exclude = seqlevels(rowRanges(counts))[grep("random|chrM|Un|GL|KI", seqlevels(rowRanges(counts)))]
if(length(chr.exclude) > 0) {
  idy = grep(paste(chr.exclude, collapse="|"), rowRanges(counts))
  if(length(idy) > 0) {counts = counts[-idy, ]}
}

# Remove peaks with 0 fragments across all cells
remove_idx <- as.numeric(which(Matrix::rowSums(assay(counts)) == 0))
if(length(remove_idx) > 0) counts <- counts[-remove_idx, ]

# ==============================================================================
# Run chromVar
# ==============================================================================

# Calculate GC bias
counts <- addGCBias(
  counts, 
  genome = BSgenome.Hsapiens.UCSC.hg38)

# Get motif annotations
anno_ix <- matchMotifs(
  pwms = pwms, 
  subject = counts, 
  genome = BSgenome.Hsapiens.UCSC.hg38)

saveRDS(anno_ix, file=paste0(out_dir, "mat.motifs_day10_COCL2_Viertra_motif_peak_match.rds"))

mat_rowSums <- rowSums(mat)
mat_colSums <- colSums(mat)

# Compute deviations
dev <- computeDeviations(
  object = counts, 
  annotations = anno_ix)

saveRDS(dev, file=paste0(out_dir, "mat.motifs_day10_COCL2_Viertra_annotations.rds"))
# All motif features
mat.motifs <- deviationScores(dev)


