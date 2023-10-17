library(tidyverse)
library(ggplot2)

load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation.RData')

sample <- 'day10_CIS'
chromVAR_correlation_with_growth <- cis_cor_vec
chromVAR_correlation_with_growth$abs_correlation <- abs(chromVAR_correlation_with_growth$correlation)
chromVAR_correlation_with_growth$p_val_adj <- p.adjust(chromVAR_correlation_with_growth$p.value, 'BH')
chromVAR_correlation_with_growth <- chromVAR_correlation_with_growth[chromVAR_correlation_with_growth$p_val_adj < 0.05, ]
chromVAR_correlation_with_growth <- chromVAR_correlation_with_growth %>%
  arrange(desc(abs_correlation))

anova_out <- read.csv(paste0('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/', sample, '_ANOVA_ChromVar_pvals.csv'))
anova_out$p_val_adj <- p.adjust(anova_out$p_val, 'BH')
anova_out_sig <- anova_out[anova_out$p_val_adj < 0.05, ]

top25 <- as.integer(nrow(chromVAR_correlation_with_growth) * 0.25)
chromVAR_correlation_with_growth_top25 <- chromVAR_correlation_with_growth[c(1: top25), ]

chromVAR_corr_growth_anova_sig <- intersect(chromVAR_correlation_with_growth_top25$motif_names,
                                            anova_out_sig$motif_names)
chromVAR_corr_growth_anova_sig <- as.data.frame(chromVAR_corr_growth_anova_sig)
colnames(chromVAR_corr_growth_anova_sig) <- c('motif_names')
chromVAR_corr_growth_anova_sig <- merge(chromVAR_corr_growth_anova_sig, chromVAR_correlation_with_growth_top25[, c('motif_names', 'correlation', 'p.value')], by='motif_names')
write.csv(chromVAR_corr_growth_anova_sig, paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/", sample, "_adaptation_motifs.csv"), row.names = FALSE)
