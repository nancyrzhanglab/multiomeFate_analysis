library(ggplot2)

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation.RData')

colnames(dabtram_cor_vec) <- c("motif_names", "correlation_dabtram", "p.value_dabtram" )
colnames(cocl2_cor_vec) <- c("motif_names", "correlation_cocl2", "p.value_cocl2" )
colnames(cis_cor_vec) <- c("motif_names", "correlation_cis", "p.value_cis" )

dabtram_cor_vec$p_val_adjust_dabtram <- p.adjust(dabtram_cor_vec$p.value_dabtram, 'BH')
cocl2_cor_vec$p_val_adjust_cocl2 <- p.adjust(cocl2_cor_vec$p.value_cocl2, 'BH')
cis_cor_vec$p_val_adjust_cis <- p.adjust(cis_cor_vec$p.value_cis, 'BH')

dabtram_cor_vec$neg_log10_pval_dabtram <- (-1) * log10(dabtram_cor_vec$p.value_dabtram)
cocl2_cor_vec$neg_log10_pval_cocl2 <- (-1) * log10(cocl2_cor_vec$p.value_cocl2)
cis_cor_vec$neg_log10_pval_cis <- (-1) * log10(cis_cor_vec$p.value_cis)

ggplot() +
  geom_point(data=dabtram_cor_vec, aes(x=correlation_dabtram, y=neg_log10_pval_dabtram), color='gray', size=1) +
  xlim(c(-0.6, 0.6)) +
  theme_bw()
ggplot() +
  geom_point(data=cocl2_cor_vec, aes(x=correlation_cocl2, y=neg_log10_pval_cocl2), color='gray', size=1) +
  xlim(c(-0.6, 0.6)) +
  theme_bw()
ggplot() +
  geom_point(data=cis_cor_vec, aes(x=correlation_cis, y=neg_log10_pval_cis), color='gray', size=1) +
  xlim(c(-0.6, 0.6)) +
  theme_bw()

hist(dabtram_cor_vec$correlation_dabtram, breaks=100)
hist(cocl2_cor_vec$correlation_cocl2, breaks=100)
hist(cis_cor_vec$correlation_cis, breaks=100)

dabtram_cor_vec <- dabtram_cor_vec[dabtram_cor_vec$p_val_adjust_dabtram < 0.05, ]
cocl2_cor_vec <- cocl2_cor_vec[cocl2_cor_vec$p_val_adjust_cocl2 < 0.05, ]
cis_cor_vec <- cis_cor_vec[cis_cor_vec$p_val_adjust_cis < 0.05, ]

corr_df <- merge(dabtram_cor_vec, cocl2_cor_vec, by='motif_names')
corr_df <- merge(corr_df, cis_cor_vec, by='motif_names')

ggplot(corr_df) +
  geom_point(aes(x = correlation_dabtram, y = correlation_cocl2), alpha=0.5) +
  xlim(c(-0.6, 0.6)) +
  ylim(c(-0.6, 0.6))

ggplot(corr_df) +
  geom_point(aes(x = correlation_dabtram, y = correlation_cis), alpha=0.5) +
  xlim(c(-0.6, 0.6)) +
  ylim(c(-0.6, 0.6))

ggplot(corr_df) +
  geom_point(aes(x = correlation_cocl2, y = correlation_cis), alpha=0.5) +
  xlim(c(-0.6, 0.6)) +
  ylim(c(-0.6, 0.6))

cor.test(corr_df$correlation_dabtram, corr_df$correlation_cocl2)
cor.test(corr_df$correlation_dabtram, corr_df$correlation_cis)
cor.test(corr_df$correlation_cocl2, corr_df$correlation_cis)
