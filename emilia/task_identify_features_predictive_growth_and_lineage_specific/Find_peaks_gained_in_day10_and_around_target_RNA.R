library(Seurat)
library(GenomicRanges)
library(stringr)
library(tidyverse)
library(ggplot2)

# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching.RData')

day0_peaks <- readRDS('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day0/Peaks/day0_peaks.rds')
day10_peaks <- readRDS('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day10_DABTRAM/Peaks/day10_DABTRAM_peaks.rds')

rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
rna_targets_non <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv')

# ==============================================================================
# Find peaks gained in day10
# ==============================================================================
day0_peaks_df <- as.data.frame(day0_peaks)
day0_peaks_df <- day0_peaks_df[, c('seqnames', 'start', 'end')]
day0_peaks <- makeGRangesFromDataFrame(day0_peaks_df)

day10_peaks_df <- as.data.frame(day10_peaks)
day10_peaks_df <- day10_peaks_df[, c('seqnames', 'start', 'end')]
day10_peaks <- makeGRangesFromDataFrame(day10_peaks_df)

overlap <- GenomicRanges::findOverlaps(day10_peaks, day0_peaks, type = 'any', select='all', ignore.strand = TRUE)
overlap_df <- as.data.frame(overlap)

# join dataframes to find chrom arms
day10_peaks_df$queryHits <- rownames(day10_peaks_df)
day0_peaks_df$subjectHits <- rownames(day0_peaks_df)

colnames(day10_peaks_df) <- c("seqnames_day10", "start_day10", "end_day10", "queryHits")

day10_peaks_df <- merge(day10_peaks_df, overlap_df, by = 'queryHits', all = TRUE)

gained_peaks <- day10_peaks_df[is.na(day10_peaks_df$subjectHits), ] 

gained_peaks <- gained_peaks[, c("seqnames_day10", "start_day10", "end_day10")]
colnames(gained_peaks) <- c('seqnames', 'start', 'end')
rownames(gained_peaks) <- seq(1: nrow(gained_peaks))

gained_peaks_GR <- makeGRangesFromDataFrame(gained_peaks)

# ==============================================================================
# Find genes that contain peaks sustained in week5
# ==============================================================================
matching_list <- Filter(function(x) length(x) > 3, matching_list)
matching_list <- Filter(function(x) length(x[["peak_names"]]) > 0, matching_list)

gene_peaks_df <- data.frame(matrix(nrow=0, ncol=2))
colnames(gene_peaks_df) <- c('gene', 'peak')
for (gene in names(matching_list)) {
  peaks <- matching_list[[gene]][["peak_names"]]
  one_gene_peak_df <- data.frame(matrix(nrow=length(peaks), ncol=2))
  colnames(one_gene_peak_df) <- c('gene', 'peak')
  one_gene_peak_df$peak <- peaks
  one_gene_peak_df$gene <- rep(gene, length(peaks))
  
  gene_peaks_df <- rbind(gene_peaks_df, one_gene_peak_df)
}


gene_peaks_df[, c('chr', 'start', 'end')] <- str_split_fixed(gene_peaks_df$peak, '-', 3)
gene_peaks_df$start <- as.integer(gene_peaks_df$start)
gene_peaks_df$end <- as.integer(gene_peaks_df$end)
gene_peaks_df <- gene_peaks_df[, c('chr', 'start', 'end', 'gene')]

gene_peaks_GR <- makeGRangesFromDataFrame(gene_peaks_df)

overlap <- GenomicRanges::findOverlaps(gene_peaks_GR, gained_peaks_GR, type = 'any', ignore.strand = TRUE)
overlap_df <- as.data.frame(overlap)

gene_peaks_df$queryHits <- rownames(gene_peaks_df)
gained_peaks$subjectHits <- rownames(gained_peaks)

colnames(gene_peaks_df) <- c("seqnames_gene", "start_gene", "end_gene", "gene", "queryHits")
colnames(gained_peaks) <- c("seqnames_gained_peaks", "start_gained_peaks", "end_gained_peaks", "subjectHits")

overlap_df <- merge(gene_peaks_df, overlap_df, by = 'queryHits')
overlap_df <- merge(gained_peaks, overlap_df, by = 'subjectHits')

overlap_df_target <- overlap_df[overlap_df$gene %in% rna_targets$gene, ]
overlap_df_nonTarget <- overlap_df[overlap_df$gene %in% rna_targets_non$gene, ]

overlap_df_target <- overlap_df_target[, c('seqnames_gained_peaks', 'start_gained_peaks', 'end_gained_peaks')]
overlap_df_nonTarget <- overlap_df_nonTarget[, c('seqnames_gained_peaks', 'start_gained_peaks', 'end_gained_peaks')]


write.table(overlap_df_target, '~/Downloads/gained_peaks_in_day10_around_target_genes.bed', sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(overlap_df_nonTarget, '~/Downloads/gained_peaks_in_day10_around_non_target_genes.bed', sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


