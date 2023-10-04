library(tidyverse)

load('~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day0/Writeup6n_lineage-imputation_day0-day10-export.RData')

sample <- 'day10_DABTRAM'

anova_out <- read.csv(paste0('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/', sample, '_processed_RNA_ANOVA_pvals.csv'))
anova_out$p_val_adj <- p.adjust(anova_out$p_val, 'BH')
anova_out_sig <- anova_out[anova_out$p_val_adj < 0.05, ]

exp_correlation_with_growth <- as.data.frame(correlation_list[[sample]])
exp_correlation_with_growth$abs_correlation <- abs(exp_correlation_with_growth$correlation)
exp_correlation_with_growth$p_val_adj <- p.adjust(exp_correlation_with_growth$p.value, 'BH')
exp_correlation_with_growth <- exp_correlation_with_growth[exp_correlation_with_growth$p_val_adj < 0.05, ]

exp_correlation_with_growth <- exp_correlation_with_growth %>%
  arrange(desc(abs_correlation))

top25 <- as.integer(nrow(exp_correlation_with_growth) * 0.25)
exp_correlation_with_growth_top25 <- exp_correlation_with_growth[c(1: top25), ]
exp_correlation_with_growth_top25$gene <- row.names(exp_correlation_with_growth_top25)

exp_corr_growth_anova_sig <- intersect(exp_correlation_with_growth_top25$gene,
                                       anova_out_sig$feature)
exp_corr_growth_anova_sig <- as.data.frame(exp_corr_growth_anova_sig)
colnames(exp_corr_growth_anova_sig) <- c('gene')
exp_corr_growth_anova_sig <- merge(exp_corr_growth_anova_sig, exp_correlation_with_growth_top25[, c('gene', 'correlation', 'p.value')], by='gene')

write.csv(exp_corr_growth_anova_sig, paste0("/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/", sample, "_adaptation_genes.csv"), row.names = FALSE)

