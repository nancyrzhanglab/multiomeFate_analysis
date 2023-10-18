library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# ==============================================================================
# Read data
# ==============================================================================
sample <- 'day10_COCL2'
# load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task5_identify_features_corr_plasticity/', sample, '_gene_exp_day10_plasticity_correlation.RData'))
load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task5_identify_features_corr_plasticity/', sample, '_motif_chromVAR_day10_plasticity_byRNA_correlation.RData'))

cor_vec$feature <- rownames(cor_vec)
cor_vec$neg_log10_pval <- (-1) * log10(cor_vec$p.value)
cor_vec$p.value_adj <- p.adjust(cor_vec$p.value, method = 'BH')
cor_vec <- cor_vec %>%
  arrange(desc(correlation))
cor_vec$order <- seq(1, nrow(cor_vec))
cor_vec <- cor_vec %>% drop_na()
# ==============================================================================
# Plotting
# ==============================================================================

features_to_label <- head(cor_vec, as.integer(nrow(cor_vec) * 0.04))$feature
features_to_label <- c(features_to_label, tail(cor_vec, as.integer(nrow(cor_vec) * 0.04))$feature)

ggplot(cor_vec, aes(x = order, y = correlation)) +
  geom_point(size=0.5) +
  # ggrepel::geom_text_repel(data = subset(cor_vec, feature %in% features_to_label),
  #                          ggplot2::aes(label = feature),
  #                          box.padding = ggplot2::unit(1, 'lines'),
  #                          point.padding = ggplot2::unit(0.1, 'lines'),
  #                          color = 'blue',
  #                          max.overlaps = 80) +
  ylab('Correlation between variance of motif chromVAR and lineage variability') +
  ylim(c(-0.8, 0.8)) +
  ggtitle(sample) +
  theme_bw()
