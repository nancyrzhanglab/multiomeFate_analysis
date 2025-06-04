rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(hdrcde)
library(ggdensity)
library(RColorBrewer)
library(circlize)
library(ggExtra)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
score_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

# =============================================================================
# Read data
# =============================================================================

fp_gini <- read.csv(paste0(output_dir, 'fate_potential_gini_index.csv'))

ada_index.DABTRAM <- read.csv(paste0(output_dir, 'adaptation_index_DABTRAM_day10_week5.csv'))
ada_index.COCL2 <- read.csv(paste0(output_dir, 'adaptation_index_COCL2_day10_week5.csv'))

ada_index.DABTRAM.minmax <- read.csv(paste0(output_dir, 'adaptation_index_DABTRAM_day10_week5_MinMax.csv'))
ada_index.COCL2.minmax <- read.csv(paste0(output_dir, 'adaptation_index_COCL2_day10_week5_MinMax.csv'))

# =============================================================================
# Wrangle
# =============================================================================
fp_gini.DABTRAM <- fp_gini[fp_gini$dataset == 'DABTRAM_d10_w5', ]
ada_index.DABTRAM <- merge(ada_index.DABTRAM, fp_gini.DABTRAM, by = 'lineage')
# ada_index.DABTRAM <- ada_index.DABTRAM[ada_index.DABTRAM$n_cells > 10, ]
ada_index.DABTRAM$adaptation.index <- (ada_index.DABTRAM$dist / mean(ada_index.DABTRAM.minmax$dist.min) ) 

fp_gini.COCL2 <- fp_gini[fp_gini$dataset == 'COCL2_d10_w5', ]
 ada_index.COCL2 <- merge(ada_index.COCL2, fp_gini.COCL2, by = 'lineage')
ada_index.COCL2$adaptation.index <- (ada_index.COCL2$dist / mean(ada_index.COCL2.minmax$dist.min))

ada_index.plot <- rbind(ada_index.DABTRAM, ada_index.COCL2)

write.csv(ada_index.plot, 
          paste0(output_dir, 'adaptation_index_norm_by_NoAdapt.csv'), 
          row.names = FALSE)

# =============================================================================
# Plot
# =============================================================================

# ggplot() +
#   # geom_point(aes(size = log10(n_cells.t2), fill = dataset), shape = 21) +
#   # geom_vline(xintercept = 0.4, linetype = 'dashed', color = 'gray') +
#   geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
#   # scale_size_continuous(range = c(1, 5)) +
#   geom_hdr(data = ada_index.DABTRAM,
#            aes(x = adaptation.index, y = gini_index, fill = after_stat(probs)),
#            color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
#   scale_fill_manual(values = brewer.pal(5, "RdPu")) +
#   geom_hdr(data = ada_index.COCL2,
#            aes(x = adaptation.index, y = gini_index, fill = after_stat(probs)),
#            color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
#   scale_fill_manual(values = brewer.pal(5, "BuGn")) +
#   # stat_cor(aes(group = dataset, color = dataset)) +
#   xlim(0, max(ada_index.plot$adaptation.index, rm.na=T)) +
#   ylim(0, 1) +
#   theme_classic() +
#   xlab('Adaptation index') +
#   ylab('Fate potential gini index')

# ada_index.plot$adaptation.index.new <- ifelse(ada_index.plot$adaptation.index.new < 0, 0, ada_index.plot$adaptation.index.new)
# ada_index.plot$adaptation.index.new <- ifelse(ada_index.plot$adaptation.index.new > 1, 1, ada_index.plot$adaptation.index.new)
p <- ggplot(ada_index.plot, aes(x = adaptation.index.new, y = gini_index, color = dataset)) +
  geom_point(aes(size = log10(n_cells.t2), fill = dataset)) +
  # geom_vline(xintercept = median(ada_index.plot$adaptation.index), linetype = 'dashed', color = 'gray') +
  # geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
  stat_cor(aes(group = dataset, color = dataset), label.x = 0.65, label.y.npc = 0.05) +
  # xlim(min(ada_index.plot$adaptation.index.new, rm.na=T), max(ada_index.plot$adaptation.index.new, rm.na=T)) +
  # xlim(0, max(ada_index.plot$adaptation.index.new, rm.na=T))+
  ylim(0, 1) +
  scale_size_continuous(range = c(1, 10), breaks=c(0, 0.5, 1, 2, 3), guide = 'none') +
  scale_fill_manual(values = c('DABTRAM_d10_w5' = '#9D85BE', 'COCL2_d10_w5' = '#6DC49C'), guide = 'none') +
  scale_color_manual(values = c('DABTRAM_d10_w5' = '#9D85BE', 'COCL2_d10_w5' = '#6DC49C'), guide = 'none') +
  theme_classic() +
  facet_wrap(~dataset, scale = 'free') +
  xlab('Adaptation index') +
  ylab('Fate potential gini index')
p <- ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
ggsave(paste0(figure_dir, 'adaptation_index_gini_index.pdf'), p, width = 6, height = 6)

p1 <- ggplot(ada_index.plot, aes(x = adaptation.index.new, y = gini_index, color = dataset)) +
  geom_point(aes(size = log10(n_cells.t2), fill = dataset)) +
  # geom_vline(xintercept = median(ada_index.plot$adaptation.index), linetype = 'dashed', color = 'gray') +
  # geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
  stat_cor(aes(group = dataset, color = dataset), label.x = 0.65, label.y.npc = 0.05) +
  xlim(min(ada_index.plot$adaptation.index.new, rm.na=T), max(ada_index.plot$adaptation.index.new, rm.na=T)) +
  ylim(0, 1) +
  scale_size_continuous(range = c(1, 10), breaks=c(0, 0.5, 1, 2, 3)) +
  scale_fill_manual(values = c('DABTRAM_d10_w5' = '#9D85BE', 'COCL2_d10_w5' = '#6DC49C')) +
  scale_color_manual(values = c('DABTRAM_d10_w5' = '#9D85BE', 'COCL2_d10_w5' = '#6DC49C')) +
  theme_classic() +
  # facet_wrap(~dataset) +
  xlab('Adaptation index') +
  ylab('Fate potential gini index')
p1
ggsave(paste0(figure_dir, 'adaptation_index_gini_index_legend.pdf'), p1, width = 6, height = 6)
