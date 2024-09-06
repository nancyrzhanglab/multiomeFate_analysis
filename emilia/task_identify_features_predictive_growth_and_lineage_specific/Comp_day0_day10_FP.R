library(ggplot2)

# ================================================================================
# Read data
# ================================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/Results_with_GP_writeup6n/day0_chromVar_day0_growth_potential_for_day10_correlation_writeup6n.RData')

# day0_cor_DABTRAM <- correlation_list[['dabtram_cor_vec']]
# day0_cor_COCL2 <- correlation_list[['cocl2_cor_vec']]
# day0_cor_CIS <- correlation_list[['cis_cor_vec']]
day0_cor_DABTRAM <- dabtram_cor_vec
day0_cor_COCL2 <- cocl2_cor_vec
day0_cor_CIS <- cis_cor_vec

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/Results_with_GP_writeup6n/day10_chromVar_day10_growth_potential_for_week5_correlation_writeup6n.RData')
# day10_cor_DABTRAM <- correlation_list[['dabtram_cor_vec']]
# day10_cor_COCL2 <- correlation_list[['cocl2_cor_vec']]
# day10_cor_CIS <- correlation_list[['cis_cor_vec']]
day10_cor_DABTRAM <- dabtram_cor_vec
day10_cor_COCL2 <- cocl2_cor_vec
day10_cor_CIS <- cis_cor_vec

# ================================================================================
# Wrangle data
# ================================================================================
# day0_cor_DABTRAM$motif <- rownames(day0_cor_DABTRAM)
# day0_cor_COCL2$motif <- rownames(day0_cor_COCL2)
# day0_cor_CIS$motif <- rownames(day0_cor_CIS)
# colnames(day0_cor_DABTRAM) <- c('correlation.day0', 'p.value.day0', 'motif')
colnames(day0_cor_DABTRAM) <- c('motif', 'correlation.day0', 'p.value.day0')
colnames(day0_cor_COCL2) <- c('motif', 'correlation.day0', 'p.value.day0')
colnames(day0_cor_CIS) <- c('motif', 'correlation.day0', 'p.value.day0')

# day10_cor_DABTRAM$motif <- rownames(day10_cor_DABTRAM)
# day10_cor_COCL2$motif <- rownames(day10_cor_COCL2)
# day10_cor_CIS$motif <- rownames(day10_cor_CIS)
# colnames(day10_cor_DABTRAM) <- c('correlation.day10', 'p.value.day10', 'motif')
colnames(day10_cor_DABTRAM) <- c('motif', 'correlation.day10', 'p.value.day10')
colnames(day10_cor_COCL2) <- c('motif', 'correlation.day10', 'p.value.day10')
colnames(day10_cor_CIS) <- c('motif', 'correlation.day10', 'p.value.day10')

# ================================================================================
# DABTRAM
# ================================================================================
comp_dabtram <- merge(day0_cor_DABTRAM, day10_cor_DABTRAM, by = 'motif')
comp_dabtram$correlation.day0 <- as.numeric(comp_dabtram$correlation.day0)
comp_dabtram$correlation.day10 <- as.numeric(comp_dabtram$correlation.day10)

tead <- comp_dabtram$motif[grepl('TEAD', comp_dabtram$motif)]
jun <- comp_dabtram$motif[grepl('JUN', comp_dabtram$motif)]
fos <- comp_dabtram$motif[grepl('FOS', comp_dabtram$motif)]
snai <- comp_dabtram$motif[grepl('SNAI', comp_dabtram$motif)]

jun_df <- comp_dabtram[comp_dabtram$motif %in% jun, ]
tead_df <- comp_dabtram[comp_dabtram$motif %in% tead, ]
snai_df <- comp_dabtram[comp_dabtram$motif %in% snai, ]

ggplot(comp_dabtram, aes(x = correlation.day0, y = correlation.day10)) +
  geom_point(alpha = 0.5) +
  stat_cor(method="pearson") +
  geom_point(data = jun_df, aes(x = correlation.day0, y = correlation.day10), color = 'red') +
  geom_point(data = tead_df, aes(x = correlation.day0, y = correlation.day10), color = 'blue') +
  geom_point(data = snai_df, aes(x = correlation.day0, y = correlation.day10), color = '#41B3A2') +
  theme_bw()

# ================================================================================
# COCL2
# ================================================================================
comp_cocl2 <- merge(day0_cor_COCL2, day10_cor_COCL2, by = 'motif')
comp_cocl2$correlation.day0 <- as.numeric(comp_cocl2$correlation.day0)
comp_cocl2$correlation.day10 <- as.numeric(comp_cocl2$correlation.day10)

tead <- comp_cocl2$motif[grepl('TEAD', comp_cocl2$motif)]
jun <- comp_cocl2$motif[grepl('JUN', comp_cocl2$motif)]
fos <- comp_cocl2$motif[grepl('FOS', comp_cocl2$motif)]
snai <- comp_cocl2$motif[grepl('SNAI', comp_cocl2$motif)]

jun_df <- comp_cocl2[comp_cocl2$motif %in% jun, ]
tead_df <- comp_cocl2[comp_cocl2$motif %in% tead, ]
snai_df <- comp_cocl2[comp_cocl2$motif %in% snai, ]

ggplot(comp_cocl2, aes(x = correlation.day0, y = correlation.day10)) +
  geom_point(alpha = 0.5) +
  stat_cor(method="pearson") +
  geom_point(data = jun_df, aes(x = correlation.day0, y = correlation.day10), color = 'red') +
  geom_point(data = tead_df, aes(x = correlation.day0, y = correlation.day10), color = 'blue') +
  geom_point(data = snai_df, aes(x = correlation.day0, y = correlation.day10), color = '#41B3A2') +
  theme_bw()

# ================================================================================
# CIS
# ================================================================================

comp_cis <- merge(day0_cor_CIS, day10_cor_CIS, by = 'motif')
comp_cis$correlation.day0 <- as.numeric(comp_cis$correlation.day0)
comp_cis$correlation.day10 <- as.numeric(comp_cis$correlation.day10)

tead <- comp_cis$motif[grepl('TEAD', comp_cis$motif)]
jun <- comp_cis$motif[grepl('JUN', comp_cis$motif)]
fos <- comp_cis$motif[grepl('FOS', comp_cis$motif)]
snai <- comp_cis$motif[grepl('SNAI', comp_cis$motif)]

jun_df <- comp_cis[comp_cis$motif %in% jun, ]
tead_df <- comp_cis[comp_cis$motif %in% tead, ]
snai_df <- comp_cis[comp_cis$motif %in% snai, ]

ggplot(comp_cis, aes(x = correlation.day0, y = correlation.day10)) +
  geom_point(alpha = 0.5) +
  stat_cor(method="pearson") +
  geom_point(data = jun_df, aes(x = correlation.day0, y = correlation.day10), color = 'red') +
  geom_point(data = tead_df, aes(x = correlation.day0, y = correlation.day10), color = 'blue') +
  geom_point(data = snai_df, aes(x = correlation.day0, y = correlation.day10), color = '#41B3A2') +
  theme_bw()
