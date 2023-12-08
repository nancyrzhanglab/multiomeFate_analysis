library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


gene_sets <- read.delim('/Users/emiliac/Dropbox/Rotations/Minn_lab/Resources/Resources genesets/Hallmark gene sets/h.all.v5.2.symbols.gmt.txt', header = FALSE)
rownames(gene_sets) <- gene_sets$V1
gene_sets <- gene_sets[ ,-which(names(gene_sets) %in% c("V1","V2"))]

EMT_genes <- unlist(gene_sets['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', ])
names(EMT_genes) <- NULL
Apoptosis_genes <- unlist(gene_sets['HALLMARK_APOPTOSIS', ])
names(Apoptosis_genes) <- NULL
Androgen_response_genes <- unlist(gene_sets['HALLMARK_ANDROGEN_RESPONSE', ])
names(Androgen_response_genes) <- NULL
Heme_genes <- unlist(gene_sets['HALLMARK_HEME_METABOLISM', ])
names(Heme_genes) <- NULL
# ==============================================================================
# Read data
# ==============================================================================
sample <- 'day10_COCL2'
# load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task5_identify_features_corr_plasticity/', sample, '_gene_exp_day10_plasticity_correlation.RData'))
# load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task5_identify_features_corr_plasticity/', sample, '_motif_chromVAR_day10_plasticity_byRNA_correlation.RData'))
load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/', sample, '_gene_exp_mean_day10_plasticity_correlation.RData'))
# load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/', sample, '_motif_mean_chromVAR_day10_plasticity_byRNA_correlation.RData'))

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

# features_to_label <- head(cor_vec, 1000)$feature
# features_to_label <- c(features_to_label, tail(cor_vec, 10)$feature)

features_to_label <- EMT_genes

ggplot(cor_vec, aes(x = order, y = correlation)) +
  geom_point(size=0.5) +
  geom_point(data = subset(cor_vec, feature %in% features_to_label), color='red') +
  # ggrepel::geom_text_repel(data = subset(cor_vec, feature %in% features_to_label),
  #                          ggplot2::aes(label = feature),
  #                          box.padding = ggplot2::unit(1, 'lines'),
  #                          point.padding = ggplot2::unit(0.1, 'lines'),
  #                          color = 'blue',
  #                          max.overlaps = 80) +
  # ylab('Corr b/w mean(chromVar) and lin. var.') +
  ylim(c(-1, 1)) +
  ggtitle(sample) +
  theme_bw()

features_to_label <- as.data.frame(features_to_label)
write.csv(features_to_label, paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/', sample, '_gene_exp_mean_day10_plasticity_correlation_TOP1000.csv'), row.names = FALSE)


cor_vec$emt_genes <- ifelse(cor_vec$feature %in% EMT_genes, 'Yes', 'No')
cor_vec$apop_genes <- ifelse(cor_vec$feature %in% Apoptosis_genes, 'Yes', 'No')
cor_vec$androgen_genes <- ifelse(cor_vec$feature %in% Androgen_response_genes, 'Yes', 'No')
cor_vec$heme_genes <- ifelse(cor_vec$feature %in% Heme_genes, 'Yes', 'No')

ggplot(cor_vec, aes(x = order, y = correlation)) +
  geom_point(data = cor_vec[cor_vec$heme_genes == 'No', ], color='gray', alpha=0.6, size=0.5) +
  geom_point(data = cor_vec[cor_vec$heme_genes == 'Yes', ], color='red', size=0.5) +
  ylab('Corr b/w mean(gene exp.) and lin. var.') +
  ylim(c(-1, 1)) +
  ggtitle(sample) +
  theme_bw()

# ==============================================================================
# Compare DABTRAM and COCL2
# ==============================================================================
# load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/day10_DABTRAM_motif_mean_chromVAR_day10_plasticity_byRNA_correlation.RData'))
# day10_dabtram <- cor_vec
# 
# load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/day10_COCL2_motif_mean_chromVAR_day10_plasticity_byRNA_correlation.RData'))
# day10_cocl2 <- cor_vec

load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/day10_DABTRAM_gene_exp_mean_day10_plasticity_correlation.RData'))
day10_dabtram <- cor_vec

load(paste0('~/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability/features_associated_with_plasticity/day10_COCL2_gene_exp_mean_day10_plasticity_correlation.RData'))
day10_cocl2 <- cor_vec


day10_dabtram$feature <- rownames(day10_dabtram)
colnames(day10_dabtram) <- c("correlation_day10_dabtram", "p.value_day10_dabtram", "feature")

day10_cocl2$feature <- rownames(day10_cocl2)
colnames(day10_cocl2) <- c("correlation_day10_cocl2", "p.value_day10_cocl2", "feature")

comp_df <- merge(day10_cocl2, day10_dabtram, by='feature')

ggplot(comp_df, aes(x = correlation_day10_dabtram, y = correlation_day10_cocl2)) +
  geom_point(alpha=0.7) + 
  ylim(c(-0.8, 0.8)) +
  xlim(c(-0.8, 0.8)) +
  theme_bw()





