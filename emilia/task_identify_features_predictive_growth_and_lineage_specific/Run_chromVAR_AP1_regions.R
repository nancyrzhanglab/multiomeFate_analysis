library(Seurat)
library(chromVAR)
library(GenomicRanges)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)

bsgenome <- BSgenome.Hsapiens.UCSC.hg38
# in_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/'
# out_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/outs/compare_chromvar_rna_corr_with_gp/'
# ==============================================================================
# Read data
# ==============================================================================
day0 <- readRDS(paste0(in_dir, "data/day0.rds"))
metadat <- day0@meta.data

ap1_peaks_targetRNA <- read.table('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_AP1_peaks.bed', sep='\t', header = TRUE)
# ap1_peaks_targetRNA <- read.table(paste0(in_dir, 'outs/compare_chromvar_rna_corr_with_gp/common_genes_in_corr_with_growth_v2_AP1_peaks.bed'), sep='\t', header = TRUE)
# ap1_peaks_targetRNA <- makeGRangesFromDataFrame(ap1_peaks_targetRNA)
ap1_peaks_targetRNA <- makeGRangesListFromDataFrame(ap1_peaks_targetRNA)
ap1_peaks_nonTargetRNA <- read.table('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget_Jun_peaks.bed', sep='\t', header = TRUE)
ap1_peaks_nonTargetRNA <- read.table(paste0(in_dir, 'outs/compare_chromvar_rna_corr_with_gp/common_genes_in_corr_with_growth_v2_nonTarget_Jun_peaks.bed'), sep='\t', header = TRUE)
ap1_peaks_nonTargetRNA <- makeGRangesFromDataFrame(ap1_peaks_nonTargetRNA)
ap1_peaks_nonTargetRNA <- makeGRangesListFromDataFrame(ap1_peaks_nonTargetRNA)

peaks <- day0@assays[["ATAC"]]@ranges

# Convert to summarized experiment
counts <- SummarizedExperiment(
  assays = list(counts = day0@assays[["ATAC"]]@counts), 
  rowRanges = peaks, 
  colData = colnames(day0))

# Remove non-standard chromosomes
chr.exclude = seqlevels(rowRanges(counts))[grep("random|chrM|Un", seqlevels(rowRanges(counts)))]
if(length(chr.exclude) > 0) {
  idy = grep(paste(chr.exclude, collapse="|"), rowRanges(counts))
  if(length(idy) > 0) {counts = counts[-idy, ]}
}

# Remove peaks with 0 fragments across all cells
remove_idx <- as.numeric(which(Matrix::rowSums(assay(counts)) == 0))
if(length(remove_idx) > 0) counts <- counts[-remove_idx, ]

# ==============================================================================
# Run chromVAR
# ==============================================================================

# Calculate GC bias
counts <- addGCBias(
  counts, 
  genome = bsgenome)

# Filter cells with insufficient reads
counts@colData@listData[["depth"]] <- metadat$nCount_ATAC
counts_filtered <- filterSamples(counts, min_depth = 1500,
                                 min_in_peaks = 0.15, shiny = FALSE)
# filtering_plot <- filterSamplesPlot(counts, min_depth = 1500, 
#                                     min_in_peaks = 0.15, use_plotly = FALSE)
# filtering_plot

# Filter peaks
counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)


# Get peak annotations
# anno_files <-makeGRangesListFromDataFrame(day0_peaks_df)
anno_files <- c(ap1_peaks_targetRNA, ap1_peaks_nonTargetRNA)
gL <- GRangesList(gr1 = ap1_peaks_targetRNA, gr2 = ap1_peaks_nonTargetRNA)
anno_names <- c('ap1_peaks_targetRNA', 'ap1_peaks_nonTargetRNA')

anno_ix <- getAnnotations(
  anno_files, 
  rowRanges=rowRanges(counts_filtered))

# Compute deviations
dev <- computeDeviations(
  object = counts_filtered, 
  annotations = anno_ix)

rownames(dev) <- anno_names
saveRDS(dev, file=paste0(outdir, "dev_day0_ap1_at_target_vs_nonTarget.rds"))