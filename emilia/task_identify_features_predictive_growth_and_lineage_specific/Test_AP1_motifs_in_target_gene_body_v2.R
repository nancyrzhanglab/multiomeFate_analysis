library(Seurat)
library(GenomicRanges)
library(dplyr)
library(tidyverse)


# ==============================================================================
# Read data
# ==============================================================================
day0 <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/day0_with_motifs.rds")
motif_info <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/motif_info.csv')
all_genes <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_all.csv')
rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
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
motif_data <- day0@assays[["ATAC"]]@motifs@data
motif_data <- motif_data[, motif_data@Dimnames[[2]] %in% tf_targets$motif_code]
motif_data <- as.data.frame(motif_data)
motif_data$num_tf_binding <- rowSums(motif_data)
motif_data <- motif_data[motif_data$num_tf_binding > 0, ]
motif_data$binding_sites <- rownames(motif_data)

binding_sites_df <- as.data.frame(motif_data$binding_sites)
colnames(binding_sites_df) <- 'binding_sites'
binding_sites_df[, c('chr', 'start', 'end')] <- str_split_fixed(binding_sites_df$binding_sites, '-', 3)

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
gene_contains_AP_motif$containAPmotif <- as.numeric(gene_contains_AP_motif$containAPmotif)
hist(gene_contains_AP_motif$containAPmotif, breaks=7, probability = TRUE)

all_genes$isTargetGene <- ifelse(all_genes$gene %in% rna_targets$gene, 'Yes', 'No')
all_genes <- merge(all_genes, gene_contains_AP_motif, by='gene')
all_genes$isContainAPMotif <- ifelse(all_genes$containAPmotif > 0, 'Yes', 'No')

summary <- all_genes %>% 
  group_by(isTargetGene, isContainAPMotif) %>% 
  summarise(n_gene = n())


all_genes$cluster_dabtram_cocl2 <- as.factor(all_genes$cluster_dabtram_cocl2)
all_genes$cluster_cis_cocl2 <- as.factor(all_genes$cluster_cis_cocl2)

summary2 <- all_genes %>% 
  group_by(cluster_dabtram_cocl2, isContainAPMotif) %>% 
  summarise(n_gene = n())
summary2_1 <- all_genes %>% 
  group_by(cluster_dabtram_cocl2) %>% 
  summarise(n = n())
summary2 <- merge(summary2, summary2_1, by = 'cluster_dabtram_cocl2')
summary2$pct <- summary2$n_gene / summary2$n * 100
summary2 <- summary2[summary2$isContainAPMotif == 'Yes', ]

all_genes_sub <- all_genes[all_genes$cluster_dabtram_cocl2 != 0 & all_genes$cluster_dabtram_cocl2 != 1 & all_genes$cluster_dabtram_cocl2 != 2, ]
ggplot(all_genes) +
  # geom_violin(aes(x = cluster_dabtram_cocl2, y = containAPmotif), scale = "width") +
  geom_jitter(aes(x = cluster_cis_cocl2, y = containAPmotif), size=0.5, alpha=0.1) +
  geom_boxplot(aes(x = cluster_cis_cocl2, y = containAPmotif), color='red', width=0.3, outlier.shape = NA)
p + stat_summary(fun.y=median, geom="point", size=2, color="red")

