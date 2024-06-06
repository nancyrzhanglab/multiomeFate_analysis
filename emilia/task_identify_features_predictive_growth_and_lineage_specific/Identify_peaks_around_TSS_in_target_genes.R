library(Seurat)
library(GenomicRanges)
library(dplyr)
library(tidyverse)

set.seed(10)
# ==============================================================================
# Read data
# ==============================================================================
day0 <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/day0_with_motifs.rds")
all_genes <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_all.csv')
rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
rna_targets_non <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv')

# tf_targets <- tf_targets[tf_targets$cluster_cis_cocl2 == 3, ]
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching_tss.RData')

# ==============================================================================
# Wrangle data
# ==============================================================================

# Get target gene body
gene_peaks <- matching_list
gene_peaks <- Filter(function(x) length(x) > 3, gene_peaks)
gene_peaks <- Filter(function(x) length(x[["peak_names"]]) > 0, gene_peaks)

gene_peaks_df <- data.frame(matrix(nrow=0, ncol=2))
colnames(gene_peaks_df) <- c('gene', 'peak')
for (gene in names(gene_peaks)) {
  peaks <- gene_peaks[[gene]][["peak_names"]]
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

gene_peaks_target <- gene_peaks_df[gene_peaks_df$gene %in% rna_targets$gene, ]
gene_peaks_non_target <- gene_peaks_df[gene_peaks_df$gene %in% rna_targets_non$gene, ]

gene_peaks_target <- gene_peaks_target[, c('chr', 'start', 'end')]
gene_peaks_non_target <- gene_peaks_non_target[, c('chr', 'start', 'end')]

write.table(gene_peaks_target, '~/Downloads/gene_peaks_tss_target.bed', sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(gene_peaks_non_target, '~/Downloads/gene_peaks_tss_non_target.bed', sep='\t',quote = FALSE, row.names = FALSE, col.names = FALSE)
