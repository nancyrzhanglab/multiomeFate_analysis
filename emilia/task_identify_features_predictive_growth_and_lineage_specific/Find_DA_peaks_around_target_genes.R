library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)

# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching.RData')

da_peaks <- read.table('~/Downloads/da_peaks_sig.bed', sep='\t')
colnames(da_peaks) <- c('chrom', 'start', 'end')

rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
rna_targets_non <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv')

# ==============================================================================
# 
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

da_peaks_GR <- makeGRangesFromDataFrame(da_peaks)

overlap <- GenomicRanges::findOverlaps(gene_peaks_GR, da_peaks_GR, type = 'any', ignore.strand = TRUE)
overlap_df <- as.data.frame(overlap)

gene_peaks_df$queryHits <- rownames(gene_peaks_df)
da_peaks$subjectHits <- rownames(da_peaks)

colnames(gene_peaks_df) <- c("seqnames_gene", "start_gene", "end_gene", "gene", "queryHits")
colnames(da_peaks) <- c("seqnames_da_peaks", "start_da_peaks", "end_da_peaks", "subjectHits")

overlap_df <- merge(gene_peaks_df, overlap_df, by = 'queryHits')
overlap_df <- merge(da_peaks, overlap_df, by = 'subjectHits')

overlap_df_target <- overlap_df[overlap_df$gene %in% rna_targets$gene, ]
overlap_df_nonTarget <- overlap_df[overlap_df$gene %in% rna_targets_non$gene, ]

overlap_df_target <- overlap_df_target[, c('seqnames_da_peaks', 'start_da_peaks', 'end_da_peaks')]
overlap_df_nonTarget <- overlap_df_nonTarget[, c('seqnames_da_peaks', 'start_da_peaks', 'end_da_peaks')]

write.table(overlap_df_target, '~/Downloads/da_peaks_in_day10_around_target_genes.bed', sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(overlap_df_nonTarget, '~/Downloads/da_peaks_in_day10_around_non_target_genes.bed', sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

