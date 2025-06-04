rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
score_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

time1 <- 'day10'
time2 <- 'week5'

# =============================================================================
# Read data
# =============================================================================

ada_index.plot <- read.csv(paste0(output_dir, 'adaptation_index_norm_by_NoAdapt.csv'))

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_COCL2.RData'))

ada_index.plot.DABTRAM <- ada_index.plot[ada_index.plot$dataset == 'DABTRAM_d10_w5', ]
ada_index.plot.COCL2 <- ada_index.plot[ada_index.plot$dataset == 'COCL2_d10_w5', ]

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.DABTRAM.time1 <- metadat[metadat$dataset == 'day10_DABTRAM', ]
metadat.DABTRAM.time2 <- metadat[metadat$dataset == 'week5_DABTRAM', ]
metadat.COCL2.time1 <- metadat[metadat$dataset == 'day10_COCL2', ]
metadat.COCL2.time2 <- metadat[metadat$dataset == 'week5_COCL2', ]

# =============================================================================
# Wrangle
# =============================================================================
umap.DABTRAM <- all_data_ft_DABTRAM_umap@cell.embeddings
umap.DABTRAM <- as.data.frame(umap.DABTRAM)
umap.DABTRAM$cell_id <- rownames(umap.DABTRAM)

umap.COCL2 <- all_data_ft_COCL2_umap@cell.embeddings
umap.COCL2 <- as.data.frame(umap.COCL2)
umap.COCL2$cell_id <- rownames(umap.COCL2)


umap.DABTRAM <- merge(umap.DABTRAM, 
                      metadat.DABTRAM.time1[, c('fatepotential_DABTRAM_d10_w5', 'cell_id')], 
                      by = 'cell_id', 
                      all.x = TRUE)


lin.DABTRAM <- 'Lin130951'  #'Lin122575'  # 'Lin27092'
adaptation.index.DABTRAM <- ada_index.plot.DABTRAM$adaptation.index[ada_index.plot.DABTRAM$lineage == lin.DABTRAM]
cells.check.time1.DABTRAM <- metadat.DABTRAM.time1[metadat.DABTRAM.time1$assigned_lineage == lin.DABTRAM, ]
cells.check.time2.DABTRAM <- metadat.DABTRAM.time2[metadat.DABTRAM.time2$assigned_lineage == lin.DABTRAM, ]

umap.cells.time1.DABTRAM <- umap.DABTRAM[umap.DABTRAM$cell_id %in% rownames(cells.check.time1.DABTRAM), ]
umap.cells.time2.DABTRAM <- umap.DABTRAM[umap.DABTRAM$cell_id %in% rownames(cells.check.time2.DABTRAM), ]

# ada_index.plot.DABTRAM <- ada_index.plot.DABTRAM[ada_index.plot.DABTRAM$n_cells.t1 >= 10, ]

ggplot(umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  # geom_point(data = umap.cells.time1, color='red', aes(fill = 'Day10')) +
  geom_point(data = umap.cells.time2.DABTRAM, color='darkgreen', aes(fill = 'Week5')) +
  geom_point(data = umap.cells.time1.DABTRAM, aes(color = fatepotential_DABTRAM_d10_w5)) +
  # ggrepel::geom_text_repel(data = umap.cells.time1, aes(label = cell_id), size = 2, max.overlaps = 20) +
  # scale_fill_manual(name="Timepoints",values = c('Day10' = 'red', 'Week5' = 'blue')) +
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(lin.DABTRAM, ' adapt. index: ', round(adaptation.index.DABTRAM, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')


example.clones <- c('Lin130951', 'Lin122575', 'Lin27092', 'Lin70082')
ggplot(ada_index.plot.DABTRAM, aes(x = log10(n_cells.t2), y = adaptation.index)) +
  geom_point(aes(size = log10(n_cells.t1)), 
              fill = '#9D85BE', shape = 21) +
  ggrepel::geom_text_repel(data = subset(ada_index.plot.DABTRAM, lineage %in% example.clones), aes(label = lineage), size = 4, max.overlaps = 20) +
  stat_cor() +
  xlim(0, 3.2) +
  ylim(0, 3) +
  ggtitle('DABTRAM') +
  scale_size_continuous(breaks=c( 1, 1.5, 2)) +
  theme_bw() +
  theme(legend.position = 'bottom' )

umap.COCL2 <- merge(umap.COCL2, 
                      metadat.COCL2.time1[, c('fatepotential_COCL2_d10_w5', 'cell_id')], 
                      by = 'cell_id', 
                      all.x = TRUE)

lin.COCL2 <- 'Lin95888'  # 'Lin39814'  #'Lin76197'  # 'Lin104509' 
adaptation.index.COCL2 <- ada_index.plot.COCL2$adaptation.index[ada_index.plot.COCL2$lineage == lin.COCL2]
cells.check.time1.COCL2 <- metadat.COCL2.time1[metadat.COCL2.time1$assigned_lineage == lin.COCL2, ]
cells.check.time2.COCL2 <- metadat.COCL2.time2[metadat.COCL2.time2$assigned_lineage == lin.COCL2, ]

umap.cells.time1.COCL2 <- umap.COCL2[umap.COCL2$cell_id %in% rownames(cells.check.time1.COCL2), ]
umap.cells.time2.COCL2 <- umap.COCL2[umap.COCL2$cell_id %in% rownames(cells.check.time2.COCL2), ]

ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  # geom_point(data = umap.cells.time1, color='red', aes(fill = 'Day10')) +
  geom_point(data = umap.cells.time2.COCL2, color='darkgreen', aes(fill = 'Week5')) +
  geom_point(data = umap.cells.time1.COCL2, aes(color = fatepotential_COCL2_d10_w5)) +
  # ggrepel::geom_text_repel(data = umap.cells.time1, aes(label = cell_id), size = 2, max.overlaps = 20) +
  # scale_fill_manual(name="Timepoints",values = c('Day10' = 'red', 'Week5' = 'blue')) +
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(lin.COCL2, ' adapt. index: ', round(adaptation.index.COCL2, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')


# ada_index.plot.COCL2 <- ada_index.plot.COCL2[ada_index.plot.COCL2$n_cells.t1 >= 10, ]

example.clones <- c('Lin95888', 'Lin39814', 'Lin76197', 'Lin104509')

ggplot(ada_index.plot.COCL2, aes(x = log10(n_cells.t2), y = adaptation.index)) +
  geom_jitter(aes(size = log10(n_cells.t1)), 
              fill = '#6DC49C', width = 0.1, height = 0.1, shape = 21) +
  stat_cor() +
  ggrepel::geom_text_repel(data = subset(ada_index.plot.COCL2, lineage %in% example.clones), aes(label = lineage), size = 4, max.overlaps = 20) +
  xlim(0, 3.5) +
  ylim(0, 1.5) +
  ggtitle('COCL2') +
  scale_size_continuous(breaks=c( 1, 1.5, 2)) +
  theme_bw() +
  theme(legend.position = 'bottom' )


# ==============================================================================
# Compare with multi-resistant clones
# ==============================================================================
metadat.week5 <- metadat[metadat$dataset %in% c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS'), ]
lin.size.week5 <- table(metadat.week5$assigned_lineage, 
                        metadat.week5$dataset)
lin.size.week5 <- as.data.frame(lin.size.week5)
lin.size.week5 <- pivot_wider(lin.size.week5, 
                              names_from = Var2, 
                              values_from = Freq)
colnames(lin.size.week5)[1] <- 'assigned_lineage'

# sum non-zero entries per row
lin.size.week5$multi_resistant <- rowSums(lin.size.week5[, c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS')] > 0)
lin.multi_resistant <- lin.size.week5[lin.size.week5$multi_resistant > 1, ] %>% pull(assigned_lineage)
ada_index.plot.DABTRAM$multi_resistant_bi <- ada_index.plot.DABTRAM$lineage %in% lin.multi_resistant

ada_index.plot.DABTRAM <- merge(ada_index.plot.DABTRAM, 
                                 lin.size.week5[, c('assigned_lineage', 'multi_resistant')], 
                                 by.x = 'lineage', by.y = 'assigned_lineage', all.x = TRUE)
ada_index.plot.DABTRAM$multi_resistant <- as.character(ada_index.plot.DABTRAM$multi_resistant)
p1.dabtram <- ggplot(ada_index.plot.DABTRAM, aes(x = log10(n_cells.t2), y = adaptation.index)) +
  geom_point(aes(size = log10(n_cells.t1), fill = multi_resistant),  shape = 21) +
  # stat_cor() +
  scale_fill_manual(values = c('1' = 'gray', '2' = 'pink', '3' = 'red')) +
  xlim(0, 3.2) +
  ylim(0, 3) +
  ggtitle('DABTRAM') +
  scale_size_continuous(breaks=c( 1, 1.5, 2)) +
  theme_bw() +
  theme(legend.position = 'bottom')

p2.dabtram <- ggplot(ada_index.plot.DABTRAM, aes(x = multi_resistant_bi, y = adaptation.index)) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape  = NA, width = 0.4) +
  geom_jitter(aes(fill = multi_resistant), shape = 21, width = 0.2, size = 2) +
  scale_fill_manual(values = c('1' = 'gray', '2' = 'pink', '3' = 'red')) +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE'))) +
  ggtitle('DABTRAM') +
  ylim(0, 3) +
  theme_bw() +
  theme(legend.position = 'bottom')

ggarrange(p1.dabtram, p2.dabtram, widths = c(1, 0.6),
          ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = 'bottom')




ada_index.plot.COCL2$multi_resistant_bi <- ada_index.plot.COCL2$lineage %in% lin.multi_resistant

ada_index.plot.COCL2 <- merge(ada_index.plot.COCL2, 
                                lin.size.week5[, c('assigned_lineage', 'multi_resistant')], 
                                by.x = 'lineage', by.y = 'assigned_lineage', all.x = TRUE)
ada_index.plot.COCL2$multi_resistant <- as.character(ada_index.plot.COCL2$multi_resistant)

p1.cocl2 <- ggplot(ada_index.plot.COCL2, aes(x = log10(n_cells.t2), y = adaptation.index)) +
  geom_point(aes(size = log10(n_cells.t1), fill = multi_resistant),  shape = 21) +
  # stat_cor() +
  scale_fill_manual(values = c('1' = 'gray', '2' = 'pink', '3' = 'red')) +
  xlim(0, 3.5) +
  ylim(0, 1.5) +
  ggtitle('COCL2') +
  scale_size_continuous(breaks=c( 1, 1.5, 2)) +
  theme_bw() +
  theme(legend.position = 'bottom')

p2.cocl2 <- ggplot(ada_index.plot.COCL2, aes(x = multi_resistant_bi, y = adaptation.index)) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape  = NA, width = 0.4) +
  geom_jitter(aes(fill = multi_resistant), shape = 21, width = 0.2, size = 2) +
  scale_fill_manual(values = c('1' = 'gray', '2' = 'pink', '3' = 'red')) +
  stat_compare_means(comparisons = list(c('FALSE', 'TRUE'))) +
  ggtitle('COCL2') +
  ylim(0, 1.5) +
  theme_bw() +
  theme(legend.position = 'bottom')

ggarrange(p1.cocl2, p2.cocl2, widths = c(1, 0.6),
          ncol = 2, nrow = 1, 
          common.legend = TRUE, legend = 'bottom')


# ==============================================================================
# Compare with adaptaiton gene expression at day10
# ==============================================================================
score_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
scores <- read.csv(paste0(score_dir, 'UCell_scores_adaptation_gene_sets.csv'), row.names = 1)
scores$dabtram.adaptation.genes_scale_UCell <- scale(scores$dabtram.adaptation.genes_UCell)
scores$cocl2.adaptation.genes_scale_UCell <- scale(scores$cocl2.adaptation.genes_UCell)

scores.DABTRAM.time1 <- scores[rownames(scores) %in% rownames(metadat.DABTRAM.time1), ]
scores.DABTRAM.time1$cell_id <- rownames(scores.DABTRAM.time1)
scores.DABTRAM.time2 <- scores[rownames(scores) %in% rownames(metadat.DABTRAM.time2), ]

metadat.DABTRAM.time1$cell_id <- rownames(metadat.DABTRAM.time1)
metadat.DABTRAM.time1 <- merge(metadat.DABTRAM.time1, scores.DABTRAM.time1, by = 'cell_id', all = TRUE)

median.scores.DABTRAM.time1 <- metadat.DABTRAM.time1 %>%
  group_by(assigned_lineage) %>%
  summarise_at(vars(ends_with('UCell')), median, na.rm = TRUE)

ada_index.plot.DABTRAM <- merge(ada_index.plot.DABTRAM, median.scores.DABTRAM.time1, by.x = 'lineage', by.y = 'assigned_lineage')
midpoint <- mean(c(scores.DABTRAM.time1$dabtram.adaptation.genes_UCell, scores.DABTRAM.time2$dabtram.adaptation.genes_UCell))


ada_index.plot.DABTRAM <- ada_index.plot.DABTRAM %>% 
  arrange(dabtram.adaptation.genes_UCell)

ggplot(ada_index.plot.DABTRAM, aes(x = adaptation.index , y = gini_index)) +
  geom_point(aes(size = log10(n_cells.t2), fill = dabtram.adaptation.genes_UCell), shape = 21) +
  scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
  scale_size_continuous(range = c(1, 8), breaks=c( 1, 2, 3)) +
  ylim(0, 1) +
  xlim(0, 2.5) +
  stat_cor() +
  theme_classic() +
  ggtitle('DABTRAM') +
  xlab('Adaptation index') +
  ylab('Fate potential gini index') +
  theme(legend.position = 'bottom')

ggplot(ada_index.plot.DABTRAM, aes(x = adaptation.index , y = dabtram.adaptation.genes_UCell)) +
  geom_point(aes(size = log10(n_cells.t2), fill = dabtram.adaptation.genes_UCell), shape = 21) +
  scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
  scale_size_continuous(range = c(1, 8), breaks=c( 1, 2, 3)) +
  stat_cor() +
  theme_classic() +
  ggtitle('DABTRAM') +
  xlab('Adaptation index') +
  ylab('Adaptation gene set UCell score') +
  theme(legend.position = 'bottom')


scores.COCL2.time1 <- scores[rownames(scores) %in% rownames(metadat.COCL2.time1), ]
scores.COCL2.time1$cell_id <- rownames(scores.COCL2.time1)
scores.COCL2.time2 <- scores[rownames(scores) %in% rownames(metadat.COCL2.time2), ]

metadat.COCL2.time1$cell_id <- rownames(metadat.COCL2.time1)
metadat.COCL2.time1 <- merge(metadat.COCL2.time1, scores.COCL2.time1, by = 'cell_id', all = TRUE)

median.scores.COCL2.time1 <- metadat.COCL2.time1 %>%
  group_by(assigned_lineage) %>%
  summarise_at(vars(ends_with('UCell')), median, na.rm = TRUE)

ada_index.plot.COCL2 <- merge(ada_index.plot.COCL2, median.scores.COCL2.time1, by.x = 'lineage', by.y = 'assigned_lineage')
midpoint <- mean(c(scores.COCL2.time1$cocl2.adaptation.genes_UCell, scores.COCL2.time2$cocl2.adaptation.genes_UCell))

ggplot(ada_index.plot.COCL2, aes(x = adaptation.index , y = gini_index)) +
  geom_point(aes(size = log10(n_cells.t2), fill = cocl2.adaptation.genes_UCell), shape = 21) +
  scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
  scale_size_continuous(range = c(1, 8), breaks=c( 1, 2, 3)) +
  ylim(0, 1) +
  xlim(0, 1.5) +
  stat_cor(label.y = 0) +
  theme_classic() +
  ggtitle('COCL2') +
  xlab('Adaptation index') +
  ylab('Fate potential gini index') +
  theme(legend.position = 'bottom')

ggplot(ada_index.plot.COCL2, aes(x = adaptation.index , y = cocl2.adaptation.genes_UCell)) +
  geom_point(aes(size = log10(n_cells.t2), fill = cocl2.adaptation.genes_UCell), shape = 21) +
  scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
  scale_size_continuous(range = c(1, 8), breaks=c( 1, 2, 3)) +
  stat_cor(method =) +
  theme_classic() +
  ggtitle('COCL2') +
  xlab('Adaptation index') +
  ylab('Adaptation gene set UCell score') +
  theme(legend.position = 'bottom')

wilcox.test(ada_index.plot.DABTRAM$adaptation.index, ada_index.plot.COCL2$adaptation.index)
wilcox.test(ada_index.plot.DABTRAM$gini_index, ada_index.plot.COCL2$gini_index)
