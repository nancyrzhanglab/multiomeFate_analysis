library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/'
# dir <- '~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
sample <- 'day10_DABTRAM'

# ==============================================================================
# Read data
# ==============================================================================
plasticity_scores <- read.table(paste0(dir, 'task0_explore_lineage_variability/outs/plasticity_scores/', sample, '_rna_d_in_v2.tsv'), sep='\t', header = TRUE)

# Read chromVars
load('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/chromVar_day10_data.RData')
# load('~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/chromVar_day10_data.RData')
chromVar_mat <- chromvar_results_dabtram

# read metadata
metadat <- read.csv(paste0(dir, sample, '/', sample, '_meta.csv'), row.names = 1)

# read motif
motifs <- read.csv(paste0(dir, 'motif_info.csv'))
# ==============================================================================
# Wrangle data
# ==============================================================================
# take only relevant columns
plasticity_scores <- plasticity_scores[, c('assigned_lineage', 'normalized_avg_euc_dist_by_shuffle')]

metadat$cell_id <- rownames(metadat)
metadat <- metadat[, c('cell_id', 'assigned_lineage')]

# Get motif scores
chromVar_mat <- as.data.frame(t(chromVar_mat))

num_motifs <- ncol(dabtram_day10_chromVar_mat)
# num_motifs <- 2

chromVar_mat$cell_id <- rownames(chromVar_mat)
chromVar_mat <- merge(chromVar_mat, metadat, by = 'cell_id')

chromVar_mat <- chromVar_mat[chromVar_mat$assigned_lineage %in% plasticity_scores$assigned_lineage, ]

# ==============================================================================
# Calculate expression variance
# ==============================================================================
chromVar_var <- data.frame(matrix(nrow=length(unique(chromVar_mat$assigned_lineage)), ncol=0))
chromVar_var$assigned_lineage <- unique(chromVar_mat$assigned_lineage)

chromVar_mat <- chromVar_mat[, -c(1)]
for (motif in colnames(chromVar_mat)[1: num_motifs]) {
  df <- chromVar_mat[, c(motif, 'assigned_lineage')]
  colnames(df) <- c('motif1', 'assigned_lineage')
  df_var <- df %>% 
    group_by(assigned_lineage) %>% 
    summarise(motif_var = var(motif1))
  colnames(df_var) <- c('assigned_lineage', motif)
  chromVar_var <- merge(chromVar_var, df_var, by = 'assigned_lineage')
}

colnames(chromVar_var) <- c('assigned_lineage', motifs$motif_names[1:2])
# ==============================================================================
# Correlate with plasticity
# ==============================================================================

chromVar_var <- merge(chromVar_var, plasticity_scores, by = 'assigned_lineage')
rownames(chromVar_var) <- chromVar_var$assigned_lineage
day10_chromVar_mat <- chromVar_var[, seq(2, ncol(chromVar_var)-1)]


cor_vec <- sapply(1:ncol(day10_chromVar_mat), function(j){
  res <- stats::cor.test(chromVar_var$normalized_avg_euc_dist_by_shuffle, day10_chromVar_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
cor_vec <- as.data.frame(t(cor_vec))
colnames(cor_vec) <- c("correlation", "p.value")
rownames(cor_vec) <- colnames(day10_chromVar_mat)

save(date_of_run, session_info,
     cor_vec,
     file = paste0(dir, "task_identify_features_predictive_growth_and_lineage_specific/outs/", sample, "_motif_chromVAR_day10_plasticity_byRNA_correlation.RData"))



