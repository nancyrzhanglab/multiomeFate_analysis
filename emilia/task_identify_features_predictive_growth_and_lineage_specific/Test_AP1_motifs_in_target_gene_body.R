library(Seurat)
library(GenomicRanges)
library(dplyr)
library(tidyverse)


# ==============================================================================
# Read data
# ==============================================================================
day0 <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/day0_with_motifs.rds")
motif_info <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/motif_info.csv')
all_genes <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/all_genes_all_corrs.csv')
rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v1.csv')
tf_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_motifs_in_corr_with_growth_v1.csv')
# load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching_tss.RData')
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching.RData')

# ==============================================================================
# Wrangle data
# ==============================================================================
# Match motif name with ID
tf_targets <- as.data.frame(tf_targets$motif_names)
colnames(tf_targets) <- 'motif_names'
tf_targets <- merge(tf_targets, motif_info, by = 'motif_names')

# Get motif binding sites
motif_data <- day0@assays[["ATAC"]]@motifs@data
motif_data <- motif_data[, motif_data@Dimnames[[2]] %in% tf_targets$motif_code]
motif_data <- as.data.frame(motif_data)
motif_data$num_tf_binding <- rowSums(motif_data)
motif_data <- motif_data[motif_data$num_tf_binding > 0, ]
motif_data$binding_sites <- rownames(motif_data)

binding_sites_df <- as.data.frame(motif_data$binding_sites)
colnames(binding_sites_df) <- 'binding_sites'
binding_sites_df[, c('chr', 'start', 'end')] <- str_split_fixed(binding_sites_df$binding_sites, '-', 3)

# Get target gene TSS
gene_peaks <- matching_list[rna_targets$gene]
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

# ==============================================================================
# Bedtool intersect (each gene separate)
# ==============================================================================
binding_sites_df <- binding_sites_df[, c("chr", "start", "end")]
binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)

gene_contains_AP_motif <- data.frame(nrow=0, ncol=2)
colnames(gene_contains_AP_motif) <- c('gene', 'containAPmotif')
for (gene in names(gene_peaks)) {
  one_gene_peaks <- gene_peaks_df[gene_peaks_df$gene == gene, ]
  one_gene_peaks <- one_gene_peaks[, c('chr', 'start', 'end')]
  one_gene_peaks_GR <- makeGRangesFromDataFrame(one_gene_peaks)
  overlap <- GenomicRanges::findOverlaps(one_gene_peaks_GR, binding_sites_GR,
                                         type = 'any', select='all', ignore.strand=TRUE)
  result <- c(gene, length(overlap))
  gene_contains_AP_motif[nrow(gene_contains_AP_motif) + 1, ] <- result
}
print(paste0('Number of genes with a peak around that contain target TF binding sites: ', nrow(gene_contains_AP_motif[gene_contains_AP_motif$containAPmotif >0, ])))


