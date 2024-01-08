library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

dir <- '~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
sample <- 'day10_DABTRAM'
genes_to_check <- c('TUBB', 'HLA-B', 'MITF')
motifs_to_check <- c('FOS::JUNB', 'MAF::NFE2', 'SOX10')
# ==============================================================================
# Read data
# ==============================================================================
plasticity_scores <- read.table(paste0(dir, sample, '/', sample, '_rna_d_in_v2.tsv'), sep='\t', header = TRUE)

# sc_obj <- readRDS(paste0(dir, 'Raw_and_Processed/', sample, '.rds'))
# rna <- sc_obj@assays[['Saver']]@data
# 
# metadat <- sc_obj@meta.data
# metadat$cell_id <- rownames(metadat)
# lineage_info <- metadat[, c('cell_id', 'assigned_lineage')]

# Read chromVars
load('~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/chromVar_day10_data.RData')
chromVar_mat <- as.data.frame(chromvar_results_dabtram)
# read metadata
metadat <- read.csv(paste0(dir, sample, '/', sample, '_meta.csv'), row.names = 1)
metadat$cell_id <- rownames(metadat)
metadat <- metadat[, c('cell_id', 'assigned_lineage')]
# read motif
motifs <- read.csv(paste0(dir, 'motif_info.csv'))

# ==============================================================================
# Wrangle data [Lineage info]
# ==============================================================================
# take only relevant columns
plasticity_scores <- plasticity_scores[, c('assigned_lineage', 'normalized_avg_euc_dist_by_shuffle')]

# calculate lineage size, filter for ones with more than 3 cells for evaluation of variance
# lineage_info_size <- plasticity_scores %>% 
#   group_by(assigned_lineage) %>% 
#   summarise(lineage_size = n())
# lineage_info_size <- lineage_info_size[lineage_info_size$lineage_size > 3, ]

# ==============================================================================
# Wrangle data [RNA]
# ==============================================================================
# Get expression
genes <- rna@Dimnames[[1]]
num_genes <- length(genes)
# num_genes <- 2

cells <- rna@Dimnames[[2]]
cells <- as.data.frame(cells)
colnames(cells) <- c('cell_id')
cells$order <- seq(1: nrow(cells))

cells <- merge(cells, lineage_info, by = 'cell_id')

rna@Dimnames[[2]] <- cells$assigned_lineage
rna <- rna[, rna@Dimnames[[2]] %in% plasticity_scores$assigned_lineage]

rna <- as.data.frame(t(as.matrix(rna)))
rna$assigned_lineage <- rownames(rna)
rna[, c('assigned_lineage', 'num')] <- str_split_fixed(rna$assigned_lineage, '\\.', 2)

# ==============================================================================
# Wrangle data [ATAC]
# ==============================================================================
chromVar_mat <- as.data.frame(t(chromVar_mat))
colnames(chromVar_mat) <- motifs$motif_names
num_motifs <- ncol(chromVar_mat)
# num_motifs <- 2

chromVar_mat$cell_id <- rownames(chromVar_mat)
chromVar_mat <- merge(chromVar_mat, metadat, by = 'cell_id')

chromVar_mat <- chromVar_mat[chromVar_mat$assigned_lineage %in% plasticity_scores$assigned_lineage, ]


# ==============================================================================
# Calculate expression variance [RNA]
# ==============================================================================
# exp_var <- data.frame(matrix(nrow=length(unique(rna$assigned_lineage)), ncol=0))
# exp_var$assigned_lineage <- unique(rna$assigned_lineage)
# 
# for (gene in genes_to_check) {
#   df <- rna[, c(gene, 'assigned_lineage')]
#   colnames(df) <- c('gene1', 'assigned_lineage')
#   df_var <- df %>% 
#     group_by(assigned_lineage) %>% 
#     summarise(exp_var = var(gene1))
#   colnames(df_var) <- c('assigned_lineage', gene)
#   exp_var <- merge(exp_var, df_var, by = 'assigned_lineage')
# }

exp_mean <- data.frame(matrix(nrow=length(unique(rna$assigned_lineage)), ncol=0))
exp_mean$assigned_lineage <- unique(rna$assigned_lineage)

for (gene in genes_to_check) {
  df <- rna[, c(gene, 'assigned_lineage')]
  colnames(df) <- c('gene1', 'assigned_lineage')
  df_mean <- df %>% 
    group_by(assigned_lineage) %>% 
    summarise(exp_mean = mean(gene1))
  colnames(df_mean) <- c('assigned_lineage', gene)
  exp_mean <- merge(exp_mean, df_mean, by = 'assigned_lineage')
}

# ==============================================================================
# Calculate expression variance [ATAC]
# ==============================================================================
# chromVar_var <- data.frame(matrix(nrow=length(unique(chromVar_mat$assigned_lineage)), ncol=0))
# chromVar_var$assigned_lineage <- unique(chromVar_mat$assigned_lineage)
# 
# chromVar_mat <- chromVar_mat[, -c(1)]
# for (motif in motifs_to_check) {
#   df <- chromVar_mat[, c(motif, 'assigned_lineage')]
#   colnames(df) <- c('motif1', 'assigned_lineage')
#   df_var <- df %>% 
#     group_by(assigned_lineage) %>% 
#     summarise(motif_var = var(motif1))
#   colnames(df_var) <- c('assigned_lineage', motif)
#   chromVar_var <- merge(chromVar_var, df_var, by = 'assigned_lineage')
# }

chromVar_mean <- data.frame(matrix(nrow=length(unique(chromVar_mat$assigned_lineage)), ncol=0))
chromVar_mean$assigned_lineage <- unique(chromVar_mat$assigned_lineage)

chromVar_mat <- chromVar_mat[, -c(1)]
for (motif in motifs_to_check) {
  df <- chromVar_mat[, c(motif, 'assigned_lineage')]
  colnames(df) <- c('motif1', 'assigned_lineage')
  df_mean <- df %>% 
    group_by(assigned_lineage) %>% 
    summarise(motif_mean = mean(motif1))
  colnames(df_mean) <- c('assigned_lineage', motif)
  chromVar_mean <- merge(chromVar_mean, df_mean, by = 'assigned_lineage')
}

# ==============================================================================
# Correlate with plasticity [RNA]
# ==============================================================================

exp_var <- merge(exp_var, plasticity_scores, by = 'assigned_lineage')

ggplot(exp_var, aes(x = normalized_avg_euc_dist_by_shuffle, y = ACTB)) +
  geom_point() +
  ylab('Var(ACTB)') +
  geom_smooth(method='lm', formula= y~x)


exp_mean <- merge(exp_mean, plasticity_scores, by = 'assigned_lineage')

colnames(exp_mean) <- c("assigned_lineage", "TUBB", "HLA_B", "MITF", "normalized_avg_euc_dist_by_shuffle")
ggplot(exp_mean, aes(x = normalized_avg_euc_dist_by_shuffle, y = HLA_B)) +
  geom_point() +
  ylab('Mean(HLA-B)') +
  geom_smooth(method='lm', formula= y~x)


# ==============================================================================
# Correlate with plasticity [ATAC]
# ==============================================================================

chromVar_var <- merge(chromVar_var, plasticity_scores, by = 'assigned_lineage')

colnames(chromVar_var) <- c("assigned_lineage", "FOS_JUN", "MAF_NFE2" , "YY1", "normalized_avg_euc_dist_by_shuffle")
ggplot(chromVar_var, aes(x = normalized_avg_euc_dist_by_shuffle, y = YY1)) +
  geom_point() +
  ylab('Var(YY1)') +
  geom_smooth(method='lm', formula= y~x)

chromVar_mean <- merge(chromVar_mean, plasticity_scores, by = 'assigned_lineage')

colnames(chromVar_mean) <- c("assigned_lineage", "FOS_JUNB", "MAF_NFE2" , "SOX10", "normalized_avg_euc_dist_by_shuffle")
ggplot(chromVar_mean, aes(x = normalized_avg_euc_dist_by_shuffle, y = SOX10)) +
  geom_point() +
  ylab('Mean(SOX10)') +
  geom_smooth(method='lm', formula= y~x)

