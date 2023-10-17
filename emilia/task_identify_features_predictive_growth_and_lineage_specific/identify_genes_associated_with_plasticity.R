library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/'
sample <- 'day10_DABTRAM'

# ==============================================================================
# Read data
# ==============================================================================
plasticity_scores <- read.table(paste0(dir, 'task0_explore_lineage_variability/outs/plasticity_scores/', sample, '_rna_d_in_v2.tsv'), sep='\t', header = TRUE)

sc_obj <- readRDS(paste0(dir, 'task_identify_features_predictive_growth_and_lineage_specific/data/', sample, '.rds'))
rna <- sc_obj@assays[['Saver']]@data

metadat <- sc_obj@meta.data
metadat$cell_id <- rownames(metadat)
lineage_info <- metadat[, c('cell_id', 'assigned_lineage')]

# ==============================================================================
# Wrangle data
# ==============================================================================
# take only relevant columns
plasticity_scores <- plasticity_scores[, c('assigned_lineage', 'normalized_avg_euc_dist_by_shuffle')]

# calculate lineage size, filter for ones with more than 3 cells for evaluation of variance
lineage_info_size <- lineage_info %>% 
  group_by(assigned_lineage) %>% 
  summarise(lineage_size = n())
lineage_info_size <- lineage_info_size[lineage_info_size$lineage_size > 3, ]

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
rna <- rna[, rna@Dimnames[[2]] %in% lineage_info_size$assigned_lineage]

rna <- as.data.frame(t(as.matrix(rna)))
rna$assigned_lineage <- rownames(rna)
rna[, c('assigned_lineage', 'num')] <- str_split_fixed(rna$assigned_lineage, '\\.', 2)

# ==============================================================================
# Calculate expression variance
# ==============================================================================
exp_var <- data.frame(matrix(nrow=length(unique(rna$assigned_lineage)), ncol=0))
exp_var$assigned_lineage <- unique(rna$assigned_lineage)

for (gene in colnames(rna)[1: num_genes]) {
  df <- rna[, c(gene, 'assigned_lineage')]
  colnames(df) <- c('gene1', 'assigned_lineage')
  df_var <- df %>% 
    group_by(assigned_lineage) %>% 
    summarise(exp_var = var(gene1))
  colnames(df_var) <- c('assigned_lineage', gene)
  exp_var <- merge(exp_var, df_var, by = 'assigned_lineage')
}

# ==============================================================================
# Correlate with plasticity
# ==============================================================================

exp_var <- merge(exp_var, plasticity_scores, by = 'assigned_lineage')
rownames(exp_var) <- exp_var$assigned_lineage
day10_exp_mat <- exp_var[, seq(2, ncol(exp_var)-1)]
  
  
cor_vec <- sapply(1:ncol(day10_exp_mat), function(j){
  res <- stats::cor.test(exp_var$normalized_avg_euc_dist_by_shuffle, day10_exp_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
cor_vec <- as.data.frame(t(cor_vec))
colnames(cor_vec) <- c("correlation", "p.value")
rownames(cor_vec) <- colnames(day10_exp_mat)

save(date_of_run, session_info,
     cor_vec,
     file = paste0(dir, "task_identify_features_predictive_growth_and_lineage_specific/outs/", sample, "_gene_exp_day10_plasticity_correlation.RData"))



