rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

res.genes <- c("WNT5A", "AXL", "EGFR", "JUN", "NGFR", "PCNA")
# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
metadat <- all_data@meta.data

# ==============================================================================
# Read correlation results (FB)
# ==============================================================================
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

load(paste0(result_dir, 'geneSaver_on_day0_cor_vec_DABTRAM.RData'))
cor_vec.DABTRAM <- cor_vec

load(paste0(result_dir, 'geneSaver_on_day0_cor_vec_COCL2.RData'))
cor_vec.COCL2 <- cor_vec

load(paste0(result_dir, 'geneSaver_on_day0_cor_vec_CIS.RData'))
cor_vec.CIS <- cor_vec

# ==============================================================================
# Read correlation results (FP)
# ==============================================================================
result_dir2 <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
load(paste0(result_dir2, 'saver_cor_vec.RData'))
dabtram_d0_saver_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 
cocl2_d0_saver_cor_vec <- saver_cor_vec[['cocl2_d0_saver_cor_vec']] 
cis_d0_saver_cor_vec <- saver_cor_vec[['cis_d0_saver_cor_vec']]

# ==============================================================================
# Wrangle
# ==============================================================================
colnames(cor_vec.DABTRAM) <- c('gene', 'cor.d0.FB', 'p_val.d0.FB')
colnames(dabtram_d0_saver_cor_vec) <- c('cor.d0.FP', 'p_val.d0.FP')
dabtram_d0_saver_cor_vec$gene <- rownames(dabtram_d0_saver_cor_vec)

df.dabtram <- merge(cor_vec.DABTRAM, dabtram_d0_saver_cor_vec, by = 'gene', all = TRUE)
df.dabtram <- df.dabtram[df.dabtram$gene %in% res.genes, ]

p1 <- ggplot(df.dabtram, aes(x = cor.d0.FB, y = cor.d0.FP)) +
  geom_point() +
  # stat_cor() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation') +
  ylab('FP correlation') +
  xlim(-0.1, 0.5) +
  ylim(-0.5, 0.1) +
  ggrepel::geom_text_repel(aes(label = gene), size = 5, box.padding = 0.5) +
  ggtitle('DABTRAM') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


colnames(cor_vec.COCL2) <- c('gene', 'cor.d0.FB', 'p_val.d0.FB')
colnames(cocl2_d0_saver_cor_vec) <- c('cor.d0.FP', 'p_val.d0.FP')
cocl2_d0_saver_cor_vec$gene <- rownames(cocl2_d0_saver_cor_vec)

df.cocl2 <- merge(cor_vec.COCL2, cocl2_d0_saver_cor_vec, by = 'gene', all = TRUE)
df.cocl2 <- df.cocl2[df.cocl2$gene %in% res.genes, ]
p2 <- ggplot(df.cocl2, aes(x = cor.d0.FB, y = cor.d0.FP)) +
  geom_point() +
  # stat_cor() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation') +
  ylab('FP correlation') +
  xlim(-0.1, 0.5) +
  ylim(-0.5, 0.1) +
  ggrepel::geom_text_repel(aes(label = gene), size = 5, box.padding = 0.5) +
  ggtitle('COCL2') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


colnames(cor_vec.CIS) <- c('gene', 'cor.d0.FB', 'p_val.d0.FB')
colnames(cis_d0_saver_cor_vec) <- c('cor.d0.FP', 'p_val.d0.FP')
cis_d0_saver_cor_vec$gene <- rownames(cis_d0_saver_cor_vec)
df.cis <- merge(cor_vec.CIS, cis_d0_saver_cor_vec, by = 'gene', all = TRUE)
df.cis <- df.cis[df.cis$gene %in% res.genes, ]
p3 <- ggplot(df.cis, aes(x = cor.d0.FB, y = cor.d0.FP)) +
  geom_point() +
  # stat_cor() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation') +
  ylab('FP correlation') +
  # xlim(-0.1, 0.5) +
  # ylim(-0.5, 0.1) +
  ggrepel::geom_text_repel(aes(label = gene), size = 5, box.padding = 0.5) +
  ggtitle('CIS') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

ggarrange(p1, p2, p3,
          ncol = 3, nrow = 1)
ggsave(paste0(figure_dir, 'Supp_correlation_d0_FatePotential_and_Features.pdf'),
       width = 9, height = 3)
