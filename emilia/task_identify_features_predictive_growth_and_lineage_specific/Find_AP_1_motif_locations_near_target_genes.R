library(Seurat)
library(GenomicRanges)
library(dplyr)
library(tidyverse)

set.seed(10)
# ==============================================================================
# Read data
# ==============================================================================
day0 <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/day0_with_motifs.rds")
motif_info <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/motif_info.csv')
all_genes <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_all.csv')
rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
rna_targets_non <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv')

tf_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_motifs_in_corr_with_growth_v2.csv')
# tf_targets <- tf_targets[tf_targets$cluster_cis_cocl2 == 3, ]
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching.RData')

# ==============================================================================
# Wrangle data
# ==============================================================================

# Match motif name with ID
tf_targets <- as.data.frame(tf_targets$motif_names)
colnames(tf_targets) <- 'motif_names'
tf_targets <- merge(tf_targets, motif_info, by = 'motif_names')

# Get motif binding sites
motif_data_all <- day0@assays[["ATAC"]]@motifs@data
motif_data <- motif_data_all[, motif_data_all@Dimnames[[2]] %in% tf_targets$motif_code]
motif_data <- as.data.frame(motif_data)
motif_data$num_tf_binding <- rowSums(motif_data)
motif_data <- motif_data[motif_data$num_tf_binding > 0, ]
motif_data$binding_sites <- rownames(motif_data)

motif_data_nonAP1 <- motif_data_all[, ! motif_data_all@Dimnames[[2]] %in% tf_targets$motif_code]
motif_data_nonAP1 <- as.data.frame(motif_data_nonAP1)
motif_data_nonAP1$num_tf_binding <- rowSums(motif_data_nonAP1)
motif_data_nonAP1 <- motif_data_nonAP1[motif_data_nonAP1$num_tf_binding > 0, ]
motif_data_nonAP1$binding_sites <- rownames(motif_data_nonAP1)

motif_data_sm <- motif_data[, c('binding_sites', 'num_tf_binding')]
colnames(motif_data_sm) <- c('binding_sites', 'num_tf_binding_AP1')
motif_data_nonAP1_sm <- motif_data_nonAP1[, c('binding_sites', 'num_tf_binding')]
colnames(motif_data_nonAP1_sm) <- c('binding_sites', 'num_tf_binding_nonAP1')
num_binding_comp <- merge(motif_data_sm, motif_data_nonAP1_sm, by = 'binding_sites', all = TRUE)
num_binding_comp[is.na(num_binding_comp)] <- 0
num_binding_comp$diff <- num_binding_comp$num_tf_binding_AP1 - num_binding_comp$num_tf_binding_nonAP1
ggplot(num_binding_comp) +
  geom_point(aes(x = num_tf_binding_AP1, y = num_tf_binding_nonAP1), size = 0.5, alpha = 0.5)

binding_sites_df <- as.data.frame(motif_data$binding_sites)
colnames(binding_sites_df) <- 'binding_sites'
binding_sites_df[, c('chr', 'start', 'end')] <- str_split_fixed(binding_sites_df$binding_sites, '-', 3)
binding_sites_df <- binding_sites_df[, c('chr', 'start', 'end')]
binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)

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

gene_peaks_GR <- makeGRangesFromDataFrame(gene_peaks_df)

# ==============================================================================
# Intersect peaks around genes and AP-1 motif locations
# ==============================================================================
overlap <- GenomicRanges::findOverlaps(binding_sites_GR, gene_peaks_GR, select="all", ignore.strand=TRUE)
overlap_df <- as.data.frame(overlap)

binding_sites_df$queryHits <- rownames(binding_sites_df) 
gene_peaks_df$subjectHits <- rownames(gene_peaks_df)

colnames(binding_sites_df) <- c("chr_TF", "start_TF", "end_TF", "queryHits")
colnames(gene_peaks_df) <- c("chr_gene", "start_gene", "end_gene", "gene", "subjectHits")

overlap_df <- merge(overlap_df, binding_sites_df, by = 'queryHits', all = TRUE)
overlap_df <- merge(overlap_df, gene_peaks_df, by = 'subjectHits', all.x = TRUE)

# ==============================================================================
# Label target versus non-target genes
# ==============================================================================
overlap_df$isTargetGene <- ifelse(overlap_df$gene %in% rna_targets$gene, 'YES', 'NO')
overlap_df$isNonTargetRNA <- ifelse(overlap_df$gene %in% rna_targets_non$gene, 'YES', 'NO')


overlap_df_targetRNA <- overlap_df[overlap_df$isTargetGene == 'YES', ] 
overlap_df_nonTargetRNA <- overlap_df[overlap_df$isNonTargetRNA == 'YES', ] 

overlap_df_targetRNA <- overlap_df_targetRNA[, c("chr_gene", "start_gene", "end_gene")]
overlap_df_nonTargetRNA <- overlap_df_nonTargetRNA[, c("chr_gene", "start_gene", "end_gene")]

write.table(overlap_df_targetRNA, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_AP1_peaks.bed', sep='\t', row.names = FALSE)
write.table(overlap_df_nonTargetRNA, '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget_AP1_peaks.bed', sep='\t', row.names = FALSE)

  
