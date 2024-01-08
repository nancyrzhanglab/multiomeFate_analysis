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
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/peak-gene-matching_tss.RData')

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
tss_sites <- matching_list[rna_targets$gene]
tss_sites <- Filter(function(x) length(x) > 3, tss_sites)
tss_sites <- lapply(tss_sites, function(x) 
  paste0(strsplit(x[["peak_names"]][1], '-', fixed=T)[[1]][1], '-',
         x[["gene_tss"]]@ranges@start))
tss_sites_df <- as.data.frame(tss_sites)
tss_sites_df <- as.data.frame(t(tss_sites_df))
colnames(tss_sites_df) <- c('TSS')
tss_sites_df$gene <- rownames(tss_sites_df)
tss_sites_df[, c('chr', 'TSS')] <- str_split_fixed(tss_sites_df$TSS, '-', 2)
tss_sites_df$TSS <- as.integer(tss_sites_df$TSS)
tss_sites_df$start <- tss_sites_df$TSS - 1000
tss_sites_df$end <- tss_sites_df$TSS + 100

# ==============================================================================
# Bedtool intersect (all together)
# ==============================================================================
binding_sites_df <- binding_sites_df[, c("chr", "start", "end")]
tss_sites_df <- tss_sites_df[, c("chr", "start", "end")]

binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)
tss_sites_GR <- makeGRangesFromDataFrame(tss_sites_df)

overlap <- GenomicRanges::findOverlaps(tss_sites_GR, binding_sites_GR,
                                       type = 'any', select='all', ignore.strand=TRUE)

# ==============================================================================
# Bedtool intersect (each motif separate)
# ==============================================================================
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
tss_contains_AP_motif <- merge(tss_contains_AP_motif, motif_info, by = 'motif_code')

# ==============================================================================
# Bedtool intersect (each gene separate)
# ==============================================================================
binding_sites_df <- binding_sites_df[, c("chr", "start", "end")]
binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)

tss_contains_AP_motif <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(tss_contains_AP_motif) <- c('gene', 'containAPmotif')

rownames(tss_sites_df) <- NULL
for (row in rownames(tss_sites_df)) {
  gene <- tss_sites_df[row, 'gene']
  one_gene <- tss_sites_df[row, c("chr", "start", "end")]
  
  one_gene_GR <- makeGRangesFromDataFrame(one_gene)
  overlap <- GenomicRanges::findOverlaps(one_gene_GR, binding_sites_GR,
                                         type = 'any', select='all', ignore.strand=TRUE)
  result <- c(gene, length(overlap))
  tss_contains_AP_motif[nrow(tss_contains_AP_motif) + 1, ] <- result
}

tss_contains_AP_motif <- merge(tss_contains_AP_motif, motif_info, by = 'motif_code')


# ==============================================================================
# Bedtool intersect (each motif separate) [Random motifs]
# ==============================================================================

# Get motif binding sites
motif_data_random <- day0@assays[["ATAC"]]@motifs@data
motif_data_random <- motif_data_random[, !(motif_data_random@Dimnames[[2]] %in% tf_targets$motif_code)]
random_index <- sample(1:motif_data_random@Dim[2], nrow(tf_targets), replace=FALSE)
motif_data_random <- motif_data_random[, random_index]
motif_data_random <- as.data.frame(motif_data_random)
motif_data_random$num_tf_binding <- rowSums(motif_data_random)
motif_data_random <- motif_data_random[motif_data_random$num_tf_binding > 0, ]
motif_data_random$binding_sites <- rownames(motif_data_random)

# bedtool intersect
tss_contains_random_motif <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(tss_contains_random_motif) <- c('motif_code', 'num_TSS_contain_binding_sites')

for (col in colnames(motif_data_random)[1: 30]) {
  one_motif <- motif_data_random[, c(col, 'binding_sites')]
  colnames(one_motif) <- c('motif', 'binding_sites')
  one_motif <- one_motif[one_motif$motif == TRUE, ]
  one_motif[, c('chr', 'start', 'end')] <- str_split_fixed(one_motif$binding_sites, '-', 3)
  one_motif <- one_motif[, c('chr', 'start', 'end')]
  one_motif_GR <- makeGRangesFromDataFrame(one_motif)
  overlap <- GenomicRanges::findOverlaps(tss_sites_GR, one_motif_GR,
                                         type = 'any', select='all', ignore.strand=TRUE)
  result <- c(col, length(overlap))
  tss_contains_random_motif[nrow(tss_contains_random_motif) + 1, ] <- result
}

tss_contains_random_motif$num_TSS_contain_binding_sites <- as.numeric(tss_contains_random_motif$num_TSS_contain_binding_sites)
tss_contains_random_motif <- merge(tss_contains_random_motif, motif_info, by = 'motif_code')


# ==============================================================================
# Fisher's exact test [Gene]
# ==============================================================================
# Get target gene TSS
tss_sites <- matching_list[all_genes$gene]
tss_sites <- Filter(function(x) length(x) > 3, tss_sites)
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
tss_sites_df$end <- tss_sites_df$TSS

tss_sites_df$isTargetRNA <- ifelse(tss_sites_df$gene %in% rna_targets$gene, 'Yes', 'No')

rownames(tss_sites_df) <- NULL
tss_sites_df$rowNum <- row.names(tss_sites_df)
tss_sites_df_sm <- tss_sites_df[, c("chr", "start", "end", "gene", "isTargetRNA")]
tss_sites_df_sm_GR <- makeGRangesFromDataFrame(tss_sites_df_sm)

binding_sites_df <- binding_sites_df[, c("chr", "start", "end")]
binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)

overlap <- GenomicRanges::findOverlaps(tss_sites_df_sm_GR, binding_sites_GR,
                                       type = 'any', select='all', ignore.strand=TRUE)

tss_sites_df$containsAPmotifs <- ifelse(tss_sites_df$rowNum %in% overlap@from, 'Yes', 'No')

contingency_table <- tss_sites_df %>% 
  group_by(isTargetRNA, containsAPmotifs) %>% 
  summarise(n_gene = n())

fisher.test(matrix(c(102,247,540,1487),nrow=2,ncol=2),alternative="greater")


# ==============================================================================
# Fisher's exact test [Motifs]
# ==============================================================================
# Get motif binding sites
motif_data <- day0@assays[["ATAC"]]@motifs@data
motif_data <- as.data.frame(motif_data)
motif_data$num_tf_binding <- rowSums(motif_data)
motif_data <- motif_data[motif_data$num_tf_binding > 0, ]
motif_data$binding_sites <- rownames(motif_data)

binding_sites_df <- as.data.frame(motif_data$binding_sites)
colnames(binding_sites_df) <- 'binding_sites'
binding_sites_df[, c('chr', 'start', 'end')] <- str_split_fixed(binding_sites_df$binding_sites, '-', 3)

# Get target gene TSS
tss_sites <- matching_list[rna_targets$gene]
tss_sites <- Filter(function(x) length(x) > 3, tss_sites)
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


binding_sites_df$isTargetTF <- ifelse(binding_sites_df$binding_sites %in% rna_targets$gene, 'Yes', 'No')

rownames(tss_sites_df) <- NULL
tss_sites_df$rowNum <- row.names(tss_sites_df)
tss_sites_df_sm <- tss_sites_df[, c("chr", "start", "end", "gene", "isTargetRNA")]
tss_sites_df_sm_GR <- makeGRangesFromDataFrame(tss_sites_df_sm)

binding_sites_df <- binding_sites_df[, c("chr", "start", "end")]
binding_sites_GR <- makeGRangesFromDataFrame(binding_sites_df)

overlap <- GenomicRanges::findOverlaps(tss_sites_df_sm_GR, binding_sites_GR,
                                       type = 'any', select='all', ignore.strand=TRUE)

tss_sites_df$containsAPmotifs <- ifelse(tss_sites_df$rowNum %in% overlap@from, 'Yes', 'No')

contingency_table <- tss_sites_df %>% 
  group_by(isTargetRNA, containsAPmotifs) %>% 
  summarise(n_gene = n())

fisher.test(matrix(c(102,247,540,1487),nrow=2,ncol=2),alternative="greater")


