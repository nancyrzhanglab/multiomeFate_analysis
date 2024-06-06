library(stringr)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)

set.seed(10)
# ==============================================================================
# Read data
# ==============================================================================
load('~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Writeup6n_correlation-with-growthpotential_day10.RData')
all_genes <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_all.csv')
rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
rna_targets_non <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv')

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching_tss.RData')

motif_data <- read.table('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/rxr_nr_d1_motifs.bed')
colnames(motif_data) <- c('Chr', 'Start', 'End', 'Name', 'Length', 'Strand')
motif_data$Start <- motif_data$Start - 50
motif_data$End <- motif_data$End + 50

# ==============================================================================
# Wrangle data
# ==============================================================================

# get correlation scores
day10_cis <- as.data.frame(cis_cor_vec)
day10_cocl2 <- as.data.frame(cocl2_cor_vec)
day10_dabtram <- as.data.frame(dabtram_cor_vec)

day10_cis$gene <- rownames(day10_cis)
day10_cocl2$gene <- rownames(day10_cocl2)
day10_dabtram$gene <- rownames(day10_dabtram)

# colnames(day10_cis) <- c("correlation.day10cis", "p.value.day10cis", "gene" )
# colnames(day10_cocl2) <- c("correlation.day10cocl2", "p.value.day10cocl2", "gene" )
# colnames(day10_dabtram) <- c("correlation.day10dabtram", "p.value.day10dabtram", "gene" )

colnames(day10_cis) <- c("correlation.day10cis", "gene" )
colnames(day10_cocl2) <- c("correlation.day10cocl2", "gene" )
colnames(day10_dabtram) <- c("correlation.day10dabtram", "gene" )

corr_df <- merge(day10_cis, day10_cocl2, by = 'gene')
corr_df <- merge(corr_df, day10_dabtram, by = 'gene')

corr_df$gene[grep('AXL', corr_df$gene)]
# corr_df_label <- corr_df[corr_df$gene %in% c("CDK1", "CDK2", "TGFB2", "FOS", "FOSL1", "JUN", "MYC","CCND1",
#                                              "EGFR", "ARF5", "MMP1", "MMP3", "CD44"), ]
corr_df_label <- corr_df[corr_df$gene %in% c("GLI2", "FOXF1", "DKK1", "WNT5A", "RUNX2", "ERBB4",
                                             "AREG", "MYC", "AXL"), ]

ggplot(corr_df, aes(x = correlation.day10cocl2, y = correlation.day10dabtram)) +
  geom_point(alpha=0.1) +
  geom_point(data = corr_df_label, color = 'red') +
  geom_text_repel(data = corr_df_label, aes(label = gene),
                  box.padding = unit(0.3, "lines"),
                  color = 'red',
                  segment.size = 0.1) +
  theme_bw()



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

# Get motif binding sites
binding_sites_df <- motif_data
binding_sites_df <- binding_sites_df[, c('Chr', 'Start', 'End')]
binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)
reduce(unlist(binding_sites_GR), min.gapwidth = 2L)

# ==============================================================================
# Intersect peaks around genes and motif locations
# ==============================================================================
overlap <- GenomicRanges::findOverlaps(binding_sites_GR, gene_peaks_GR, select="all", ignore.strand=TRUE)
overlap_df <- as.data.frame(overlap)

binding_sites_df$queryHits <- rownames(binding_sites_df) 
gene_peaks_df$subjectHits <- rownames(gene_peaks_df)

colnames(binding_sites_df) <- c("chr_TF", "start_TF", "end_TF", "queryHits")
colnames(gene_peaks_df) <- c("chr_gene", "start_gene", "end_gene", "gene", "subjectHits")

overlap_df <- merge(overlap_df, binding_sites_df, by = 'queryHits')
overlap_df <- merge(overlap_df, gene_peaks_df, by = 'subjectHits')

# ==============================================================================
# Label target versus non-target genes
# ==============================================================================
overlap_df$isTargetGene <- ifelse(overlap_df$gene %in% rna_targets$gene, 'YES', 'NO')
overlap_df$isNonTargetRNA <- ifelse(overlap_df$gene %in% rna_targets_non$gene, 'YES', 'NO')

# overlap_df$label <- ifelse(overlap_df$gene %in% rna_targets$gene, 'in', 'NO')
# ==============================================================================
# Label target versus non-target genes
# ==============================================================================
corr_df$hasRXR <- ifelse(corr_df$gene %in% overlap_df$gene, 'Contains RXR(NR)_D1 motif', 'Does not contain RXR(NR)_D1')
ggplot(corr_df) +
  geom_point(aes(x = correlation.day10cocl2, y = correlation.day10dabtram, col = hasRXR), alpha=0.5) +
  facet_wrap(. ~ hasRXR) +
  scale_color_manual(values=c('#f084b8', '#5db8f5')) +
  theme_bw()
