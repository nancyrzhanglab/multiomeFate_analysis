library(tidyverse)
library(data.table)

dir <- '~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
future_treatment <- 'day10_DABTRAM'

# ==============================================================================
# Read data
# ==============================================================================

# growth potential
load(paste0(dir, future_treatment, '/Writeup6n_DABTRAM_day10_lineage-imputation_stepdown_concise-postprocessed.RData'))
growth_potential_use <- cell_imputed_count
growth_potential_use <- growth_potential_use[!is.na(growth_potential_use)]

# adaptation motifs
adaptation_motifs <- read.csv(paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/Results_with_GP_writeup6n/", future_treatment, "_adaptation_motifs.csv"))

# day0 data
file_name <- 'chromVar_day0_data'
load(paste0(dir, 'ChromVAR/JASPAR/', file_name, '.RData'))
data <- chromvar_results_day0
data <- as.data.frame(data)

metadata_day0 <- read.csv(paste0(dir, 'day0/day0_meta.csv'), row.names = 1)
metadata_day0$cell_id <- rownames(metadata_day0)

# day10 data
metadata_day10 <- read.csv(paste0(dir, future_treatment, '/', future_treatment, '_meta.csv'), row.names = 1)
metadata_day10$cell_id <- rownames(metadata_day10)

motifs <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/motif_info.csv')

data$motif_code <- rownames(as.data.frame(data))
data <- merge(data, motifs, by='motif_code')
rownames(data) <- data$motif_names
data <- select(data, -c('motif_code', 'motif_names'))
# ==============================================================================
# Wrangle data
# ==============================================================================

# =====================================
# find winning lineages (v1 - max(day10 growth potential for week5) > 0)
# =====================================
# winners <- as.data.frame(growth_potential_use[growth_potential_use > 0])
# others <- as.data.frame(growth_potential_use[growth_potential_use <= 0])
# 
# colnames(winners) <- c('growth_potential')
# winners$category  <- 'winners'
# winners$cell_id <- row.names(winners)
# winners <- merge(winners, metadata_day10[, c('assigned_lineage', 'cell_id')], by='cell_id')
# 
# metadata_day10$winning_lineage <- ifelse(metadata_day10$assigned_lineage %in% winners$assigned_lineage, 'winning_lineage', 'other_lineage')
# 
# day10_lineage_status <- metadata_day10[, c('assigned_lineage', 'winning_lineage')]
# row.names(day10_lineage_status) <- NULL
# day10_lineage_status <- unique( day10_lineage_status )
# 
# day10_winners <- day10_lineage_status[day10_lineage_status$winning_lineage == 'winning_lineage', ]
# day10_others <- day10_lineage_status[day10_lineage_status$winning_lineage == 'other_lineage', ]

# =====================================
# find winning lineages (v2 - max(day10 growth potential for week5) > 0 + mean(day10 growth potential for week5) > -1)
# =====================================

growth_potential_use <- as.data.frame(growth_potential_use)
colnames(growth_potential_use) <- c('growth_potential')
growth_potential_use$cell_id <- row.names(growth_potential_use)
growth_potential_use <- merge(growth_potential_use, metadata_day10[, c('assigned_lineage', 'cell_id')], by='cell_id')

day10_lineage_status <- growth_potential_use %>% 
  group_by(assigned_lineage) %>% 
  summarise(max_GP = max(growth_potential),
            mean_GP = mean(growth_potential),
            lineage_size = n())

day10_lineage_status$winning_lineage <- ifelse(day10_lineage_status$max_GP > 0 & day10_lineage_status$mean_GP > -1,
                                          'winning_lineage', 'other_lineage')

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
# looking at raw data
# =====================================
data_winners_all <- as.data.frame(data[, winners_df$cell_id])
data_others_all <- as.data.frame(data[, others_df$cell_id])

data_winners_all <- as.data.frame(t(data_winners_all))
data_others_all <- as.data.frame(t(data_others_all))

data_winners_all$category <- 'winners'
data_others_all$category <- 'others'

data_all <- as.data.frame(rbind(data_winners_all, data_others_all))

adaptation_motifs <- adaptation_motifs[adaptation_motifs$motif_names %in% colnames(data_all), ]
adaptation_motifs$abs_correlation <- abs(adaptation_motifs$correlation)
adaptation_motifs$corr_dir <- ifelse(adaptation_motifs$correlation > 0, 'Pos', 'None')
adaptation_motifs$corr_dir <- ifelse(adaptation_motifs$correlation < 0, 'Neg', adaptation_motifs$corr_dir)
adaptation_motifs <- adaptation_motifs %>%
  arrange(desc(abs_correlation))
top_adaptation_motifs <- head(adaptation_motifs, 20)
bottom_adaptation_motifs <- tail(adaptation_motifs, 5)
other_motif_names <- adaptation_motifs %>%
  filter(! motif_names %in% top_adaptation_motifs$motif_names) %>%
  filter(! motif_names %in% bottom_adaptation_motifs$motif_names) %>%
  sample_n(10)

top_adaptation_motifs$motif_names_category <- 'Top'
bottom_adaptation_motifs$motif_names_category <- 'Bottom'
other_motif_names$motif_names_category <- 'Others'

# motif_names_to_check <- rbind(top_adaptation_motifs, bottom_adaptation_motifs, other_motif_names)
# motif_names_to_check <- top_adaptation_motifs
motif_names_to_check <- adaptation_motifs

data_to_check <- data_all[, c(motif_names_to_check$motif_names, 'category')]
data_to_check$cell_id <- row.names(data_to_check)
data_to_check <- melt(data_to_check, id_vars=c('cell_id', 'category'))
colnames(data_to_check) <- c('cell_lineage_category', 'cell_id', 'motif_names', 'chromVAR')

data_to_check <- merge(data_to_check, motif_names_to_check, by='motif_names')


# =====================================
# calculate the difference in expression level between winner and other cells
# =====================================
data_mean <- merge(data_winners, data_others, by='feature')
data_mean$winner_minus_others <- data_mean$winners_mean - data_mean$others_mean
data_mean <- merge(data_mean, motifs, by.x = 'feature', by.y = 'motif_names')

results <- merge(adaptation_motifs, data_mean, by.x = 'motif_names', by.y = 'feature', all=TRUE)
results$corr_dir <- ifelse(results$correlation > 0, 'Pos', 'None')
results$corr_dir <- ifelse(results$correlation < 0, 'Neg', results$corr_dir)

results[is.na(results)] <- 1
results$corr_dir <- ifelse(results$corr_dir == 1, 'Others', results$corr_dir)
results <- results[results$winner_minus_others != 1, ]
# results <- results[results$correlation != 1, ]

# =====================================
# Performing t-test
# =====================================
columns <- c('motif_names', 'mean_winner', 'mean_other', 't_statistic', 'p_val')
t_test_results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(t_test_results) <- columns

for (m in motif_names_to_check$motif_names) {
  m_winner <- data_to_check[data_to_check$motif_names == m & data_to_check$cell_lineage_category == 'winners', ]
  m_other <- data_to_check[data_to_check$motif_names == m & data_to_check$cell_lineage_category == 'others', ]
  m_t_test <- t.test(m_winner$chromVAR, m_other$chromVAR, alternative = 'two.sided')
  t_statistics <- m_t_test[["statistic"]][["t"]]
  t_test_p_val <- m_t_test[["p.value"]] 
  one_motif <- data.frame(matrix(ncol = 5, nrow = 1))
  colnames(one_motif) <- columns
  one_motif$t_statistic <- t_statistics
  one_motif$p_val <- t_test_p_val
  one_motif$mean_winner <- mean(m_winner$chromVAR)
  one_motif$mean_other <- mean(m_other$chromVAR)
  one_motif$motif_names <- m
  
  t_test_results <- rbind(t_test_results, one_motif)
}

t_test_results$p_val_adjust <- p.adjust(t_test_results$p_val, 'BH')
t_test_results$neg_log10_pval <- (-1) * log10(t_test_results$p_val)

# =====================================
# annotate 
# =====================================
t_test_results$winner_minus_others <- t_test_results$mean_winner - t_test_results$mean_other

# ==============================================================================
#  Plotting
# ==============================================================================

ggplot(results, aes(x = factor(corr_dir, levels = c('Pos', 'Neg', 'Others')), y=winner_minus_others)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(alpha=0.3, size=1) +
  xlab('Correlation direction') +
  ylim(c(-0.7, 0.7)) +
  theme_bw()

results$JUN_related <- ifelse(grepl('JUN', results$motif_names), 'Yes', 'No')

ggplot(results, aes(x = factor(corr_dir, levels = c('Pos', 'Neg', 'Others')), y=winner_minus_others)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(alpha=0.3, size=1, aes(color = JUN_related)) +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Correlation direction') +
  ggtitle('Day0 for DABTRAM') +
  ylim(c(-0.7, 0.7)) +
  theme_bw()



ggplot(data_to_check, aes(x = reorder(motif_names, -abs_correlation), y = chromVAR, fill=cell_lineage_category, color=corr_dir)) +
  geom_boxplot(lwd=0.8) +
  scale_fill_manual(values=c("#A7C7E7", "#ff0000")) +
  scale_color_manual(values=c("#A9A9A9", "#50C878")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust = 0.2, hjust=0.2)) +
  coord_cartesian(ylim = c(-5, 5))

# =====================================
# plot volcano plots
# =====================================
t_test_results <- merge(t_test_results, adaptation_motifs[, c("motif_names", "correlation")], by='motif_names')
t_test_results <- t_test_results %>%
  arrange(desc(correlation))
features_to_label <- head(t_test_results, 10)$motif_names
features_to_label <- c(features_to_label, tail(t_test_results, 10)$motif_names)

sig_threshold <- min(t_test_results[t_test_results$p_val_adjust < 0.05, ]$neg_log10_pval)
ggplot(data=t_test_results, aes(x=winner_minus_others, y=neg_log10_pval)) +
  geom_hline(yintercept = sig_threshold, linetype='dashed', alpha=0.7) +
  geom_point(color='gray', size=1) +
  ggrepel::geom_text_repel(data = subset(t_test_results, motif_names %in% features_to_label),
                           ggplot2::aes(label = motif_names),
                           box.padding = ggplot2::unit(1, 'lines'),
                           point.padding = ggplot2::unit(0.1, 'lines'),
                           color = 'blue',
                           max.overlaps = 80) +
  xlim(c(-1.7, 1.7)) +
  theme_minimal()

hist(t_test_results$p_val, breaks=100)
