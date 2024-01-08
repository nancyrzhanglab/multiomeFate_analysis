library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)


# ==============================================================================
# Read data
# ==============================================================================
# load('~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day0/Writeup6n_lineage-imputation_day0-day10-export.RData')
# load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day10_chromVar_day10_growth_potential_for_week5_correlation.RData')
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/day0_chromVar_day0_growth_potential_for_day10_correlation.RData')

# modality <-'Gene_exp'
modality <- 'ChromVar'
# sample <- 'day10_COCL2'
sample <- 'day0_CIS'
# anova_out <- read.csv(paste0('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/', sample, '_processed_RNA_ANOVA_pvals.csv'))
# anova_out <- read.csv(paste0('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day0_processed_RNA_ANOVA_pvals.csv'))
# anova_out <- read.csv(paste0('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/', sample, '_ANOVA_', modality, '_pvals.csv'))
anova_out <- read.csv(paste0('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task1_day10_week5_lineage_specific/day0_ANOVA_', modality, '_pvals.csv'))
anova_out$feature <- anova_out$motif_names

anova_out$p_val_adj_anova <- p.adjust(anova_out$p_val, 'BH')
anova_out$neg_log10_p_val <- (-1) * log10(anova_out$p_val)

# exp_correlation_with_growth <- as.data.frame(correlation_list[[sample]])
exp_correlation_with_growth <- cis_cor_vec
colnames(exp_correlation_with_growth) <- c( "feature", "correlation","p.value")
exp_correlation_with_growth$abs_correlation <- abs(exp_correlation_with_growth$correlation)
exp_correlation_with_growth$p.value_adj_corr <- p.adjust(exp_correlation_with_growth$p.value, 'BH')


# ==============================================================================
# Pick significant features
# ==============================================================================
exp_correlation_with_growth$feature <- rownames(exp_correlation_with_growth)

to_plot <- merge(anova_out, exp_correlation_with_growth, by = 'feature')
to_plot$anova_sig <- ifelse(to_plot$p_val_adj_anova < 0.05, 'Yes', 'No')
to_plot$corr_sig <- ifelse(to_plot$p.value_adj_corr < 0.05, 'Yes', 'No')

to_plot$neg_log10_p_val <- ifelse(to_plot$neg_log10_p_val > 200, 200, to_plot$neg_log10_p_val)
anova_threshold <- min(to_plot[to_plot$p_val_adj_anova < 0.05, ]$neg_log10_p_val)

to_plot <- to_plot %>%
  arrange(desc(correlation))

num_features <- nrow(to_plot)
num_lineage_specific_adaptation_genes <- nrow(to_plot[to_plot$anova_sig == 'Yes' & 
                                                      to_plot$corr_sig == 'Yes', ])

title <- paste0(sample, ' (', num_lineage_specific_adaptation_genes, ' number of sig-corr. + lineage-spec. TF binding, out of ', num_features, ' )')
genes_of_interest <- head(to_plot, as.integer(nrow(to_plot) * 0.01))$feature
genes_of_interest <- c(genes_of_interest, tail(to_plot, as.integer(nrow(to_plot) * 0.01))$feature)
ggplot(to_plot, aes(x = correlation, y = neg_log10_p_val, color = corr_sig)) +
  geom_point(size = 1, alpha=0.6) +
  ggrepel::geom_text_repel(data = subset(to_plot, feature %in% genes_of_interest),
                           ggplot2::aes(label = feature),
                           box.padding = ggplot2::unit(1, 'lines'),
                           point.padding = ggplot2::unit(0.1, 'lines'),
                           color = 'blue',
                           max.overlaps = 80) +
  scale_color_manual(values=c('black','red')) +
  geom_hline(yintercept = anova_threshold, linetype = "dashed", color='#696969') +
  ylab('neg_log_10_p_val_ANOVA') +
  xlim(c(-0.8, 0.8)) +
  ggtitle(title) +
  theme_bw()
ggsave(paste0('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/', modality, '_corr_w_growth_potential_vs_anova_p_', sample, '.png'),
dpi = 300)
