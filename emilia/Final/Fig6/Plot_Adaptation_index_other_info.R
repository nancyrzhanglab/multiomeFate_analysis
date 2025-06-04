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
# Plot DABTRAM
# =============================================================================
umap.DABTRAM <- all_data_ft_DABTRAM_umap@cell.embeddings
umap.DABTRAM <- as.data.frame(umap.DABTRAM)
umap.DABTRAM$cell_id <- rownames(umap.DABTRAM)

umap.COCL2 <- all_data_ft_COCL2_umap@cell.embeddings
umap.COCL2 <- as.data.frame(umap.COCL2)
umap.COCL2$cell_id <- rownames(umap.COCL2)


umap.DABTRAM <- merge(umap.DABTRAM, 
                      metadat.DABTRAM.time1[, c('fatepotential_DABTRAM_d10_w5', 'cell_id', 'assigned_lineage')], 
                      by = 'cell_id', 
                      all.x = TRUE)


example.clones <- c('Lin130951', 'Lin122575')
adaptation.index.DABTRAM.lin1 <- ada_index.plot.DABTRAM$adaptation.index[ada_index.plot.DABTRAM$lineage == example.clones[1]]
cells.check.time1.DABTRAM.lin1 <- metadat.DABTRAM.time1[metadat.DABTRAM.time1$assigned_lineage == example.clones[1], ]
cells.check.time2.DABTRAM.lin1 <- metadat.DABTRAM.time2[metadat.DABTRAM.time2$assigned_lineage == example.clones[1], ]

umap.cells.time1.DABTRAM.lin1 <- umap.DABTRAM[umap.DABTRAM$cell_id %in% rownames(cells.check.time1.DABTRAM.lin1), ]
umap.cells.time2.DABTRAM.lin1 <- umap.DABTRAM[umap.DABTRAM$cell_id %in% rownames(cells.check.time2.DABTRAM.lin1), ]


p1.dabtram <- ggplot(umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  geom_point(size = 0, aes(fill = fatepotential_DABTRAM_d10_w5), shape = NA) +
  geom_point(data = umap.cells.time2.DABTRAM.lin1, color='#623594',  shape = 18) +
  geom_point(data = umap.cells.time1.DABTRAM.lin1, aes(fill = fatepotential_DABTRAM_d10_w5), shape =21, size = 2) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(example.clones[1], ' adapt. index: ', round(adaptation.index.DABTRAM.lin1, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')

adaptation.index.DABTRAM.lin2 <- ada_index.plot.DABTRAM$adaptation.index[ada_index.plot.DABTRAM$lineage == example.clones[2]]
cells.check.time1.DABTRAM.lin2 <- metadat.DABTRAM.time1[metadat.DABTRAM.time1$assigned_lineage == example.clones[2], ]
cells.check.time2.DABTRAM.lin2 <- metadat.DABTRAM.time2[metadat.DABTRAM.time2$assigned_lineage == example.clones[2], ]

umap.cells.time1.DABTRAM.lin2 <- umap.DABTRAM[umap.DABTRAM$cell_id %in% rownames(cells.check.time1.DABTRAM.lin2), ]
umap.cells.time2.DABTRAM.lin2 <- umap.DABTRAM[umap.DABTRAM$cell_id %in% rownames(cells.check.time2.DABTRAM.lin2), ]


p2.dabtram <- ggplot(umap.DABTRAM, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  geom_point(size = 0, aes(fill = fatepotential_DABTRAM_d10_w5), shape = NA) +
  geom_point(data = umap.cells.time2.DABTRAM.lin2, color='#623594',  shape = 18) +
  geom_point(data = umap.cells.time1.DABTRAM.lin2, aes(fill = fatepotential_DABTRAM_d10_w5), shape =21, size = 2) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(example.clones[2], ' adapt. index: ', round(adaptation.index.DABTRAM.lin2, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')

p3.dabtram <- ggplot(ada_index.plot.DABTRAM, aes(x = log10(n_cells.t2), y = adaptation.index)) +
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

p.dabtram <- ggarrange(p3.dabtram, 
          ggarrange(p2.dabtram, p1.dabtram, ncol = 1, nrow = 2, common.legend = TRUE, legend = 'bottom'),
          ncol = 2, nrow = 1, common.legend = TRUE, legend = 'bottom', widths = c(1, 0.6))


ggsave(p.dabtram, 
         filename = paste0(figure_dir, 'Supp_DABTRAM_adaptation_index.pdf'), 
         width = 8, height = 6)


# =============================================================================
# Plot CoCl2
# =============================================================================

umap.COCL2 <- merge(umap.COCL2, 
                    metadat.COCL2.time1[, c('fatepotential_COCL2_d10_w5', 'cell_id')], 
                    by = 'cell_id', 
                    all.x = TRUE)
example.clones <- c('Lin39814', 'Lin104509')

adaptation.index.COCL2.lin1 <- ada_index.plot.COCL2$adaptation.index[ada_index.plot.COCL2$lineage == example.clones[1]]
cells.check.time1.COCL2.lin1 <- metadat.COCL2.time1[metadat.COCL2.time1$assigned_lineage == example.clones[1], ]
cells.check.time2.COCL2.lin1 <- metadat.COCL2.time2[metadat.COCL2.time2$assigned_lineage == example.clones[1], ]

umap.cells.time1.COCL2.lin1 <- umap.COCL2[umap.COCL2$cell_id %in% rownames(cells.check.time1.COCL2.lin1), ]
umap.cells.time2.COCL2.lin1 <- umap.COCL2[umap.COCL2$cell_id %in% rownames(cells.check.time2.COCL2.lin1), ]

p1.cocl2 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  geom_point(size = 0, aes(fill = fatepotential_COCL2_d10_w5), shape = NA) +
  geom_point(data = umap.cells.time2.COCL2.lin1, color='#0F8241',  shape = 18) +
  geom_point(data = umap.cells.time1.COCL2.lin1, aes(fill = fatepotential_COCL2_d10_w5), shape =21, size = 2) +
  # ggrepel::geom_text_repel(data = umap.cells.time1, aes(label = cell_id), size = 2, max.overlaps = 20) +
  # scale_fill_manual(name="Timepoints",values = c('Day10' = 'red', 'Week5' = 'blue')) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(example.clones[1], ' adapt. index: ', round(adaptation.index.COCL2.lin1, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')

adaptation.index.COCL2.lin2 <- ada_index.plot.COCL2$adaptation.index[ada_index.plot.COCL2$lineage == example.clones[2]]
cells.check.time1.COCL2.lin2 <- metadat.COCL2.time1[metadat.COCL2.time1$assigned_lineage == example.clones[2], ]
cells.check.time2.COCL2.lin2 <- metadat.COCL2.time2[metadat.COCL2.time2$assigned_lineage == example.clones[2], ]

umap.cells.time1.COCL2.lin2 <- umap.COCL2[umap.COCL2$cell_id %in% rownames(cells.check.time1.COCL2.lin2), ]
umap.cells.time2.COCL2.lin2 <- umap.COCL2[umap.COCL2$cell_id %in% rownames(cells.check.time2.COCL2.lin2), ]

p2.cocl2 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  geom_point(size = 0, aes(fill = fatepotential_COCL2_d10_w5), shape = NA) +
  geom_point(data = umap.cells.time2.COCL2.lin2, color='#0F8241',  shape = 18) +
  geom_point(data = umap.cells.time1.COCL2.lin2, aes(fill = fatepotential_COCL2_d10_w5), shape =21, size = 2) +
  # ggrepel::geom_text_repel(data = umap.cells.time1, aes(label = cell_id), size = 2, max.overlaps = 20) +
  # scale_fill_manual(name="Timepoints",values = c('Day10' = 'red', 'Week5' = 'blue')) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(example.clones[1], ' adapt. index: ', round(adaptation.index.COCL2.lin2, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')

p3.cocl2 <- ggplot(ada_index.plot.COCL2, aes(x = log10(n_cells.t2), y = adaptation.index)) +
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


p.cocl2 <- ggarrange(p3.cocl2, 
                       ggarrange(p2.cocl2, p1.cocl2, ncol = 1, nrow = 2, common.legend = TRUE, legend = 'bottom'),
                       ncol = 2, nrow = 1, common.legend = TRUE, legend = 'bottom', widths = c(1, 0.6))


ggsave(p.cocl2, 
       filename = paste0(figure_dir, 'Supp_COCL2_adaptation_index.pdf'), 
       width = 8, height = 6)

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
# 
# ggplot(ada_index.plot.DABTRAM, aes(x = adaptation.index , y = gini_index)) +
#   geom_point(aes(size = log10(n_cells.t2), fill = dabtram.adaptation.genes_UCell), shape = 21) +
#   scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
#   scale_size_continuous(range = c(1, 8), breaks=c( 1, 2, 3)) +
#   ylim(0, 1) +
#   xlim(0, 2.5) +
#   stat_cor() +
#   theme_classic() +
#   ggtitle('DABTRAM') +
#   xlab('Adaptation index') +
#   ylab('Fate potential gini index') +
#   theme(legend.position = 'bottom')

p.adaptGEX <- ggplot(ada_index.plot.DABTRAM, aes(x = adaptation.index , y = dabtram.adaptation.genes_UCell)) +
  geom_point(aes(size = log10(n_cells.t2), fill = dabtram.adaptation.genes_UCell), shape = 21) +
  scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
  scale_size_continuous(range = c(1, 8), breaks=c( 1, 2, 3)) +
  stat_cor() +
  theme_classic() +
  ggtitle('DABTRAM') +
  xlab('Adaptation index') +
  ylab('Adaptation gene set UCell score') +
  theme(legend.position = 'bottom')

ggsave(p.adaptGEX, 
       filename = paste0(figure_dir, 'Supp_DABTRAM_adaptation_gex.pdf'), 
       width = 4.5, height = 5)
