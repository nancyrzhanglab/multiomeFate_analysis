rm(list=ls())
library(dplyr)
library(ggplot2)

# =============================================================================
# Load data
# =============================================================================

TIME_CUR = 'day10' # 'day10', or 'day0'
TIME_FUT = 'week5' # 'week5', 'day10'
TREATMENT = 'CIS' # 'COCL2', 'DABTRAM', or 'CIS'
MODALITY = 'saver_sample' # or 'peakvi' 
SAMPLE_NAME = paste0(TIME_CUR, '_', TREATMENT)

data_dir = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/"

lineage_variability = read.csv(paste0(data_dir, SAMPLE_NAME, '/', 'lineage_variability_', SAMPLE_NAME, '_saver_sample.csv'))
load('~/nzhanglab/project/Multiome_fate/out/kevin/Writeup10a/Writeup10a_data_empty.RData')
metadat = all_data@meta.data

# =============================================================================
# Wrangle data
# =============================================================================
# calculate lineage size
lineage_barcodes = metadat %>% 
  group_by(assigned_lineage, dataset) %>% 
  summarise(n_cells = n()) %>% 
  arrange(desc(n_cells))

lineage_barcodes_FUT = lineage_barcodes[lineage_barcodes$dataset == paste0(TIME_FUT, '_', TREATMENT), ]
colnames(lineage_barcodes_FUT) = c('assigned_lineage', 'dataset', 'n_cells_FUT')

lineage_variability = merge(lineage_variability, lineage_barcodes_FUT, by = 'assigned_lineage', all.x = TRUE)
lineage_variability$n_cells_FUT_LOG10 = log10(lineage_variability$n_cells_FUT + 1)
# =============================================================================
# Plotting
# =============================================================================

cor <- cor.test(lineage_variability$normalized_avg_eud_dist_by_shuffle, lineage_variability$n_cells_FUT_LOG10, method = 'spearman')
p_val = cor$p.value
rho = cor$estimate

p = ggplot(lineage_variability, aes(x = normalized_avg_eud_dist_by_shuffle, y = n_cells_FUT_LOG10)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  labs(y = 'Lineage size at future timepoint (log10)', x = 'Lineage variability') + 
  ggtitle(paste0(SAMPLE_NAME, '\n', 'rho = ', round(rho, 2), ', p = ', p_val)) +
  theme_bw()

ggsave(paste0(data_dir, SAMPLE_NAME, '/Lineage_variability_vs_future_lineage_size_', SAMPLE_NAME, '.png'), plot = p, width = 5, height = 5 )
