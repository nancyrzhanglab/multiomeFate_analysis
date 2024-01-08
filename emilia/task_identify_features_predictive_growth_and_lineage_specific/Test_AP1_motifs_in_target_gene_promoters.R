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
rna_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
tf_targets <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_motifs_in_corr_with_growth_v2.csv')
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching_tss.RData')
# load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching.RData')

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
tss_sites <- matching_list
tss_sites <- Filter(function(x) length(x) > 3, tss_sites)
tss_sites <- Filter(function(x) length(x[["peak_names"]]) > 0, tss_sites)
tss_sites <- lapply(tss_sites, function(x) 
                      paste0(strsplit(x[["peak_names"]][1], '-', fixed=T)[[1]][1], '-',
                             x[["gene_tss"]]@ranges@start))
tss_sites_df <- as.data.frame(tss_sites)
tss_sites_df <- as.data.frame(t(tss_sites_df))
colnames(tss_sites_df) <- c('TSS')
tss_sites_df$gene <- rownames(tss_sites_df)
tss_sites_df[, c('chr', 'TSS')] <- str_split_fixed(tss_sites_df$TSS, '-', 2)
tss_sites_df$TSS <- as.integer(tss_sites_df$TSS)
tss_sites_df$start <- tss_sites_df$TSS - 500
tss_sites_df$end <- tss_sites_df$TSS + 100
tss_sites_df$isTargetGene <- ifelse(tss_sites_df$gene %in% rna_targets$gene, 'Yes', 'No')

# ==============================================================================
# Bedtool intersect (each gene separate)
# ==============================================================================
binding_sites_df <- binding_sites_df[, c("chr", "start", "end")]
binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)

tss_sites_random_df <- tss_sites_df[tss_sites_df$isTargetGene == 'No', ]

tss_sites_random_sample_df <- sample_n(tss_sites_random_df, nrow(rna_targets))
tss_sites_random_sample_df <- tss_sites_random_sample_df[, c("chr", "start", "end", "gene", "TSS", "isTargetGene")]

tss_random_contains_AP_motif <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(tss_random_contains_AP_motif) <- c('gene', 'containAPmotif')

tss_all_contains_AP_motif <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(tss_all_contains_AP_motif) <- c('gene', 'containAPmotif')


rownames(tss_sites_random_sample_df) <- NULL
rownames(tss_sites_df) <- NULL
for (row in rownames(tss_sites_df)) {
  gene <- tss_sites_df[row, 'gene']
  one_gene <- tss_sites_df[row, c("chr", "start", "end")]
  
  one_gene_GR <- makeGRangesFromDataFrame(one_gene)
  overlap <- GenomicRanges::findOverlaps(one_gene_GR, binding_sites_GR,
                                         type = 'any', select='all', ignore.strand=TRUE)
  result <- c(gene, length(overlap))
  tss_all_contains_AP_motif[nrow(tss_all_contains_AP_motif) + 1, ] <- result
}
print(paste0('Number of genes whose promoter contain target TF binding sites: ', nrow(tss_all_contains_AP_motif[tss_all_contains_AP_motif$containAPmotif >0, ])))
print(paste0('Total number of genes: ', nrow(tss_all_contains_AP_motif)))

tss_all_contains_AP_motif$isTargetGene <- ifelse(tss_all_contains_AP_motif$gene %in% rna_targets$gene, 'Yes', 'No')
tss_all_contains_AP_motif$isContainAP1motif <- ifelse(tss_all_contains_AP_motif$containAPmotif > 0 , 'Yes', 'No')

tss_all_contains_AP_motif <- tss_all_contains_AP_motif[tss_all_contains_AP_motif$gene %in% all_genes$gene, ]
summary <- tss_all_contains_AP_motif %>% 
  group_by(isTargetGene, isContainAP1motif) %>% 
  summarise(n = n())
dat <- data.frame(
  "has_AP" = c(113, 620),
  "no_AP" = c(207, 1327),
  row.names = c("target_gene", "Non-target_gene"),
  stringsAsFactors = FALSE
)
chisq.test(dat)

tss_sites_df <- tss_sites_df[, c("chr", "start", "end")]
tss_sites_GR <- makeGRangesFromDataFrame(tss_sites_df)

tss_contains_AP_motif <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(tss_contains_AP_motif) <- c('motif_code', 'num_TSS_contain_binding_sites')

for (col in colnames(motif_data)[1: 27]) {
  one_motif <- motif_data[, c(col, 'binding_sites')]
  colnames(one_motif) <- c('motif', 'binding_sites')
  one_motif <- one_motif[one_motif$motif == TRUE, ]
  one_motif[, c('chr', 'start', 'end')] <- str_split_fixed(one_motif$binding_sites, '-', 3)
  one_motif <- one_motif[, c('chr', 'start', 'end')]
  one_motif_GR <- makeGRangesFromDataFrame(one_motif)
  overlap <- GenomicRanges::findOverlaps(tss_sites_GR, one_motif_GR,
                                         type = 'any', select='all', ignore.strand=TRUE)
  result <- c(col, length(overlap))
  tss_contains_AP_motif[nrow(tss_contains_AP_motif) + 1, ] <- result
}

tss_contains_AP_motif$num_TSS_contain_binding_sites <- as.numeric(tss_contains_AP_motif$num_TSS_contain_binding_sites)
