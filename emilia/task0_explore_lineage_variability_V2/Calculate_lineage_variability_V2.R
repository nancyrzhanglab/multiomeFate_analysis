rm(list=ls())
library(multiomeFate)
library(dplyr)
library(ggplot2)
# library(ggpubr)

source('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/git/multiomeFate/R/data_loader.R')
# =============================================================================
# Load data
# =============================================================================
# args = commandArgs(trailingOnly=TRUE)
# TIME = args[1]
# TREATMENT = args[2]

TIME = 'day10' # 'week5', 'day10', or 'day0'
TREATMENT = 'COCL2' # 'COCL2', 'DABTRAM', or 'CIS'
MODALITY = 'peakvi' # or 'saver_treatment'

if (TIME == 'day10') {
  FP_NAME = paste0('fatepotential_', TREATMENT, '_d10_w5')
  SAMPLE_NAME = paste0(TIME, '_', TREATMENT)
}else if (TIME == 'day0') {
  FP_NAME = paste0('fatepotential_', TREATMENT, '_d0_d10')
  SAMPLE_NAME = TIME
}else if (TIME == 'week5') {
  SAMPLE_NAME = paste0(TIME, '_', TREATMENT)
}else {
  stop('TIME not recognized')
}


# data_dir = '~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
output_dir = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/"

all_data = multiomeFate:::data_loader(which_files = c("lineage", MODALITY))
metadat = all_data@meta.data

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# =============================================================================
# Wrangle data
# =============================================================================

if (MODALITY == 'saver_treatment') {
  data_use = all_data@reductions[[paste0('Saver.', TREATMENT, '.pca')]]
} else if (MODALITY == 'peakvi') {
  data_use = all_data@reductions[[paste0('peakVI.', TREATMENT)]]
} else {
  stop('MODALITY not recognized')
}

data_use = as.data.frame(data_use@cell.embeddings)

NUM_DIM = dim(data_use)[2]
MIN_CELL = 4

# subset to only the cells in the dataset
metadat_use = metadat[metadat$dataset == SAMPLE_NAME, ]
metadat_use$barcode = rownames(metadat_use)

data_use = data_use[rownames(metadat_use), ]

data_use$barcode = rownames(data_use)
data_use = merge(data_use, metadat_use, by = 'barcode')

# calculate lineage size
lineage_barcodes = data_use %>% 
  group_by(assigned_lineage) %>% 
  summarise(n_cells = n()) %>% 
  arrange(desc(n_cells))


# =============================================================================
# Calculate pairwise distance within a lineage barcode 
# =============================================================================

d_in_df = data.frame(matrix(nrow = 0, ncol = 3))
colnames(d_in_df) = c('assigned_lineage', 'avg_euc_dist', 'n_cells')

if (MODALITY == 'saver_treatment') {
  embedding_col_names = paste0('Saver', TREATMENT, 'PC_', 1:NUM_DIM)
} else if (MODALITY == 'peakvi') {
  embedding_col_names = paste0('peakVI', TREATMENT, '_', 1:NUM_DIM)
} else {
  stop('MODALITY not recognized')
}

for(l in unique(data_use$assigned_lineage)) {
  df = data_use %>% filter(assigned_lineage == l)
  if(nrow(df) >= MIN_CELL) {
    dists = dist(df[, embedding_col_names], method = "euclidean")
    avg_euc_dist = mean(dists)
  }else {
    avg_euc_dist = 0
  }
  d_in_df = rbind(d_in_df, 
                  data.frame(assigned_lineage = l, 
                             avg_euc_dist = avg_euc_dist, 
                             n_cells = nrow(df)))
  
}

d_in_df = d_in_df[d_in_df$n_cells >= MIN_CELL, ]

# =============================================================================
# Calculate pairwise distance shuffled barcodes
# =============================================================================
data_use$shuffled_assigned_lineage = sample(data_use$assigned_lineage, size = nrow(data_use), replace = F)

d_in_shuffled_df = data.frame(matrix(nrow = 0, ncol = 3))
colnames(d_in_shuffled_df) = c('assigned_lineage', 'avg_euc_dist', 'n_cells')

for(l in unique(data_use$assigned_lineage)) {
  df = data_use %>% filter(shuffled_assigned_lineage == l)
  if(nrow(df) >= MIN_CELL) {
    dists = dist(df[, embedding_col_names], method = "euclidean")
    avg_euc_dist = mean(dists)
  }else {
    avg_euc_dist = 0
  }
  d_in_shuffled_df = rbind(d_in_shuffled_df, 
                           data.frame(assigned_lineage = l, 
                                      avg_euc_dist = avg_euc_dist, 
                                      n_cells = nrow(df)))
}

d_in_shuffled_df = d_in_shuffled_df[d_in_shuffled_df$n_cells >= MIN_CELL, ]

# =============================================================================
# Process Results
# =============================================================================
print(dim(d_in_df))
print(dim(d_in_shuffled_df))

d_in_df$normalized_avg_eud_dist_by_shuffle = d_in_df$avg_euc_dist / median(d_in_shuffled_df$avg_euc_dist)
d_in_shuffled_df$normalized_avg_eud_dist_by_shuffle = d_in_shuffled_df$avg_euc_dist / median(d_in_shuffled_df$avg_euc_dist)

write.csv(d_in_df, paste0(output_dir, SAMPLE_NAME, '/lineage_variability_', SAMPLE_NAME, '_', MODALITY, '.csv'), row.names = F)
write.csv(d_in_shuffled_df, paste0(output_dir, SAMPLE_NAME, '/lineage_variability_shuffled', SAMPLE_NAME,'_', MODALITY, '.csv'), row.names = F)

# =============================================================================
# Plotting
# =============================================================================

# Plotting lineage size distribution
p1 = ggplot(lineage_barcodes, aes(x = n_cells)) + 
  geom_histogram(binwidth = 5) + 
  theme_bw() +
  labs(title = paste0('Lineage size distribution (', SAMPLE_NAME, ', ', MODALITY, ')'), y = 'Number of cells')

ggsave(paste0(output_dir, SAMPLE_NAME, '/lineage_size_histogram_', SAMPLE_NAME, '_', MODALITY, '.png'), plot = p1, width = 5, height = 4 )

# Plotting lineage variability in histogram
d_in_df$category = 'Data'
d_in_shuffled_df$category = 'Shuffled'

d_in_plot = rbind(d_in_df, d_in_shuffled_df)
p2 = ggplot(d_in_plot, aes(x = normalized_avg_eud_dist_by_shuffle)) + 
  geom_histogram(aes(y = ..density.., fill = category), col = 'black', alpha = 0.5, position = 'identity') +
  geom_density(aes(col = category), lwd = 1) +
  theme_bw() +
  labs(title = paste0('Lineage variability (', SAMPLE_NAME, ', ', MODALITY, ')'), 
       x = 'Normalized average Euclidean distance')
ggsave(paste0(output_dir, SAMPLE_NAME, '/normalized_avg_eud_dist_', SAMPLE_NAME, '_', MODALITY, '.png'), plot = p2, width = 6, height = 4 )

# Plotting lineage variability's correlation with lineage size
p3 = ggplot(d_in_df, aes(x = n_cells, y = normalized_avg_eud_dist_by_shuffle)) + 
  geom_point() + 
  theme_bw() +
  labs(title = paste0('Lineage variability vs lineage size (', SAMPLE_NAME, ', ', MODALITY, ')'), 
       x = 'Number of cells', y = 'Normalized average Euclidean distance')
ggsave(paste0(output_dir, SAMPLE_NAME, '/lineage_variability_vs_size_', SAMPLE_NAME, '_', MODALITY, '.png'), plot = p3, width = 5, height = 5 )

# Plotting lineage variability's correlation with FP variance

if (TIME == 'day10') {
  fp_var = metadat_use %>%
    group_by(assigned_lineage) %>%
    summarise(fp_var = var(.data[[FP_NAME]])) %>%
    arrange(desc(fp_var))
  
  d_in_df = merge(d_in_df, fp_var, by = 'assigned_lineage')
  
  p4 = ggplot(d_in_df, aes(x = fp_var, y = normalized_avg_eud_dist_by_shuffle)) +
    geom_point() +
    # stat_cor(method="spearman") +
    theme_bw() +
    labs(title = paste0('Lineage variability vs FP variance (', SAMPLE_NAME, ', ', MODALITY, ')'),
         x = 'FP variance', y = 'Normalized average Euclidean distance')
  ggsave(paste0(output_dir,SAMPLE_NAME, '/lineage_variability_vs_fp_var_', SAMPLE_NAME, '_', MODALITY, '.png'),
         plot = p4, width = 5, height = 5)
  
}
