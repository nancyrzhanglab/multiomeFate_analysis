rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggplot2)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'


remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

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

dabtram_d10_saver_cor_vec <- saver_cor_vec[['dabtram_d10_saver_cor_vec']] 
cocl2_d10_saver_cor_vec <- saver_cor_vec[['cocl2_d10_saver_cor_vec']] 
cis_d10_saver_cor_vec <- saver_cor_vec[['cis_d10_saver_cor_vec']]

# ==============================================================================
# Wrangle
# ==============================================================================
colnames(cor_vec.DABTRAM) <- c('gene', 'cor.d0.FB', 'p_val.d0.FB')
colnames(dabtram_d0_saver_cor_vec) <- c('cor.d0.FP', 'p_val.d0.FP')
dabtram_d0_saver_cor_vec$gene <- rownames(dabtram_d0_saver_cor_vec)

df <- merge(cor_vec.DABTRAM, dabtram_d0_saver_cor_vec, by = 'gene', all = TRUE)
ggplot(df, aes(x = cor.d0.FB, y = cor.d0.FP)) +
  geom_point() +
  # stat_cor() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation') +
  ylab('FP correlation') +
  ggtitle('DABTRAM') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

colnames(dabtram_d10_saver_cor_vec) <- c('cor.d10.FP', 'p_val.d10.FP')
dabtram_d10_saver_cor_vec$gene <- rownames(dabtram_d10_saver_cor_vec)
df <- merge(cor_vec.DABTRAM, dabtram_d10_saver_cor_vec, by = 'gene', all = TRUE)
df <- df[order(df$cor.d0.FB, decreasing = T), ]
df$cor.d0.FB.rank <- seq(1, nrow(df), 1)
df <- df[order(df$cor.d10.FP, decreasing = T), ]
df$cor.d10.FP.rank <- seq(1, nrow(df), 1)

ggplot(df, aes(x = cor.d0.FB.rank, y = cor.d10.FP.rank)) +
  geom_point() +
  # stat_cor() +
  # geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  # geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation rank') +
  ylab('FP correlation rank') +
  ggtitle('DABTRAM') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# dbscan cluster on df
library(dbscan)
dbscan_res <- dbscan(df[, c('cor.d0.FB.rank', 'cor.d10.FP.rank')], eps = 170, minPts = 50)
df$cluster <- dbscan_res$cluster
ggplot(df, aes(x = cor.d0.FB.rank, y = cor.d10.FP.rank, color = as.factor(cluster))) +
  geom_point() +
  # stat_cor() +
  # geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  # geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation rank') +
  ylab('FP correlation rank') +
  ggtitle('DABTRAM') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

df.0 <- df[df$cluster == '0', ]

colnames(cor_vec.COCL2) <- c('gene', 'cor.d0.FB', 'p_val.d0.FB')
colnames(cocl2_d0_saver_cor_vec) <- c('cor.d0.FP', 'p_val.d0.FP')
cocl2_d0_saver_cor_vec$gene <- rownames(cocl2_d0_saver_cor_vec)
df <- merge(cor_vec.COCL2, cocl2_d0_saver_cor_vec, by = 'gene', all = TRUE)
ggplot(df, aes(x = cor.d0.FB, y = cor.d0.FP)) +
  geom_point() +
  # stat_cor() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation') +
  ylab('FP correlation') +
  ggtitle('COCL2') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

colnames(cocl2_d10_saver_cor_vec) <- c('cor.d10.FP', 'p_val.d10.FP')
cocl2_d10_saver_cor_vec$gene <- rownames(cocl2_d10_saver_cor_vec)

df <- merge(cor_vec.COCL2, cocl2_d10_saver_cor_vec, by = 'gene', all = TRUE)
df <- df[order(df$cor.d0.FB, decreasing = T), ]
df$cor.d0.FB.rank <- seq(1, nrow(df), 1)
df <- df[order(df$cor.d10.FP, decreasing = T), ]
df$cor.d10.FP.rank <- seq(1, nrow(df), 1)

ggplot(df, aes(x = cor.d0.FB, y = cor.d10.FP)) +
  geom_point() +
  # stat_cor() +
  # geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  # geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation rank') +
  ylab('FP correlation rank') +
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
df <- merge(cor_vec.CIS, cis_d0_saver_cor_vec, by = 'gene', all = TRUE)
ggplot(df, aes(x = cor.d0.FB, y = cor.d0.FP)) +
  geom_point() +
  # stat_cor() +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation') +
  ylab('FP correlation') +
  ggtitle('CIS') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

colnames(cis_d10_saver_cor_vec) <- c('cor.d10.FP', 'p_val.d10.FP')
cis_d10_saver_cor_vec$gene <- rownames(cis_d10_saver_cor_vec)
df <- merge(cor_vec.CIS, cis_d10_saver_cor_vec, by = 'gene', all = TRUE)
df <- df[order(df$cor.d0.FB, decreasing = T), ]
df$cor.d0.FB.rank <- seq(1, nrow(df), 1)
df <- df[order(df$cor.d10.FP, decreasing = T), ]
df$cor.d10.FP.rank <- seq(1, nrow(df), 1)
ggplot(df, aes(x = cor.d0.FB.rank, y = cor.d10.FP.rank)) +
  geom_point() +
  # stat_cor() +
  # geom_vline(xintercept = 0, linetype = 'dashed', color = 'red') +
  # geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  xlab('FB correlation rank') +
  ylab('FP correlation rank') +
  ggtitle('CIS') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
