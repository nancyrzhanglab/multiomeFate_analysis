library(tidyverse)
library(data.table)

dir <- '~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
future_treatment <- 'day10_COCL2'

# ==============================================================================
# Read data
# ==============================================================================

# growth potential
load(paste0(dir, future_treatment, '/Writeup6n_COCL2_day10_lineage-imputation_stepdown_concise-postprocessed.RData'))
growth_potential_use <- cell_imputed_count
growth_potential_use <- growth_potential_use[!is.na(growth_potential_use)]

# adaptation motifs
adaptation_motifs <- read.csv(paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/", future_treatment, "_adaptation_motifs.csv"))

# day0 data
file_name <- 'chromVar_day0_data'
load(paste0(dir, file_name, '.RData'))
data <- chromvar_results_day0

metadata_day0 <- read.csv(paste0(dir, 'day0/day0_meta.csv'), row.names = 1)
metadata_day0$cell_id <- rownames(metadata_day0)

# day10 data
metadata_day10 <- read.csv(paste0(dir, future_treatment, '/', future_treatment, '_meta.csv'), row.names = 1)
metadata_day10$cell_id <- rownames(metadata_day10)

motifs <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/motif_info.csv')
# ==============================================================================
# Wrangle data
# ==============================================================================

# =====================================
# find winning lineages
# =====================================
winners <- as.data.frame(growth_potential_use[growth_potential_use > 0])
others <- as.data.frame(growth_potential_use[growth_potential_use <= 0])

colnames(winners) <- c('growth_potential')
winners$category  <- 'winners'
winners$cell_id <- row.names(winners)
winners <- merge(winners, metadata_day10[, c('assigned_lineage', 'cell_id')], by='cell_id')

metadata_day10$winning_lineage <- ifelse(metadata_day10$assigned_lineage %in% winners$assigned_lineage, 'winning_lineage', 'other_lineage')

day10_lineage_status <- metadata_day10[, c('assigned_lineage', 'winning_lineage')]
row.names(day10_lineage_status) <- NULL
day10_lineage_status <- unique( day10_lineage_status )

day10_winners <- day10_lineage_status[day10_lineage_status$winning_lineage == 'winning_lineage', ]
day10_others <- day10_lineage_status[day10_lineage_status$winning_lineage == 'other_lineage', ]

# =====================================
# pick out cells at day0 in winner vs other lineages 
# =====================================
metadata_day0 <- metadata_day0[metadata_day0$assigned_lineage %in% day10_lineage_status$assigned_lineage, ]

metadata_day0$day10_lineage_status <- ifelse(metadata_day0$assigned_lineage %in% day10_winners$assigned_lineage, 'day10_winning_lineage', 'NA')
metadata_day0$day10_lineage_status <- ifelse(metadata_day0$assigned_lineage %in% day10_others$assigned_lineage, 'day10_other_lineage', metadata_day0$day10_lineage_status)

# metadata_day0_day10_lineage_status_summary <- metadata_day0 %>%
#   group_by(day10_lineage_status) %>%
#   summarise(size = n_distinct(cell_id))


winners_df <- metadata_day0[metadata_day0$day10_lineage_status == 'day10_winning_lineage', ]
others_df <- metadata_day0[metadata_day0$day10_lineage_status == 'day10_other_lineage', ]

data_winners <- as.data.frame(data[, winners_df$cell_id])
data_others <- as.data.frame(data[, others_df$cell_id])

# =====================================
# calculate the mean expression levels
# =====================================
data_winners$winners_mean <- rowMeans(data_winners, dims = 1)
data_others$others_mean <- rowMeans(data_others, dims = 1)

data_winners$feature <- rownames(data_winners)
data_others$feature <- rownames(data_others)

data_winners <- data_winners[, c('feature', 'winners_mean')]
data_others <- data_others[, c('feature','others_mean')]

# =====================================
# calculate the difference in expression level between winner and other cells
# =====================================
data_mean <- merge(data_winners, data_others, by='feature')
data_mean$winner_minus_others <- data_mean$winners_mean - data_mean$others_mean
data_mean <- merge(data_mean, motifs, by.x = 'feature', by.y = 'motif_code')

results <- merge(adaptation_motifs, data_mean, by = 'motif_names', all=TRUE)
results$corr_dir <- ifelse(results$correlation > 0, 'Pos', 'None')
results$corr_dir <- ifelse(results$correlation < 0, 'Neg', results$corr_dir)

results[is.na(results)] <- 1
results$corr_dir <- ifelse(results$corr_dir == 1, 'Others', results$corr_dir)
results <- results[results$winner_minus_others != 1, ]
# results <- results[results$correlation != 1, ]

# ==============================================================================
#  Plotting
# ==============================================================================

ggplot(results, aes(x = factor(corr_dir, levels = c('Pos', 'Neg', 'Others')), y=winner_minus_others)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(alpha=0.3, size=1) +
  xlab('Correlation direction') +
  ylim(c(-0.7, 0.7)) +
  theme_bw()

