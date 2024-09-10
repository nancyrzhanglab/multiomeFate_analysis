rm(list=ls())
library(dplyr)
library(ggplot2)

# =============================================================================
# Load data
# =============================================================================

TIME = 'week5' # 'week5', 'day10', or 'day0'
TREATMENT = 'COCL2' # 'COCL2', 'DABTRAM', or 'CIS'
MODALITY_1 = 'saver_sample' 
MODALITY_2 = 'peakvi'
SAMPLE_NAME = paste0(TIME, '_', TREATMENT)

data_dir = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/"


lineage_variability_m1 = read.csv(paste0(data_dir, SAMPLE_NAME, '/', 'lineage_variability_', SAMPLE_NAME, '_saver_sample.csv'))
lineage_variability_m2 = read.csv(paste0(data_dir, SAMPLE_NAME, '/embeding_treatment/', 'lineage_variability_', SAMPLE_NAME, '_', MODALITY_2, '.csv'))

stopifnot(dim(lineage_variability_m1)[1] == dim(lineage_variability_m2)[1])

# =============================================================================
# Wrangle data
# =============================================================================
colnames(lineage_variability_m1) = c('assigned_lineage', paste0('avg_euc_dist_', MODALITY_1), 
                                     'n_cells', paste0('normalized_avg_eud_dist_by_shuffle_', MODALITY_1))

colnames(lineage_variability_m2) = c('assigned_lineage', paste0('avg_euc_dist_', MODALITY_2), 
                                     'n_cells', paste0('normalized_avg_eud_dist_by_shuffle_', MODALITY_2))

to_plot = merge(lineage_variability_m1, lineage_variability_m2, by = c('assigned_lineage', 'n_cells'))

# plot correlation between lineage variability by both modalities
# to_plot = to_plot[to_plot$n_cells >= 10, ]
cor = cor.test(to_plot[, paste0('normalized_avg_eud_dist_by_shuffle_', MODALITY_1)], 
               to_plot[, paste0('normalized_avg_eud_dist_by_shuffle_', MODALITY_2)], method = 'pearson')
p_val = cor$p.value
rho = cor$estimate

p = ggplot(to_plot, aes(x = get(paste0('normalized_avg_eud_dist_by_shuffle_', MODALITY_1)), 
                        y = get(paste0('normalized_avg_eud_dist_by_shuffle_', MODALITY_2)))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(x = paste0('Average Euclidean distance by ', MODALITY_1), y = paste0('Average Euclidean distance by ', MODALITY_2)) +
  ggtitle(paste0(SAMPLE_NAME, '\n', 'rho = ', round(rho, 2), ', p = ', p_val)) +
  theme_bw()

ggsave(paste0(data_dir, SAMPLE_NAME, '/RNA_ATAC_Comp', SAMPLE_NAME, '.png'), plot = p, width = 5, height = 5 )



