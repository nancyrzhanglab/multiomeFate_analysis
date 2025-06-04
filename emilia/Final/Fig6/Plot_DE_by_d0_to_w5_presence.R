rm(list=ls())

library(tidyverse)
library(ggplot2)

results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

keygenes <- c("WNT5A", "AXL", "EGFR", "JUN", "NGFR", "PCNA")
# ==============================================================================
# Read data
# ==============================================================================

de.dabtram <- read.csv(paste0(results_dir, 'day0_DE_genes_by_w5DABTRAM_presence.csv'))
de.cocl2 <- read.csv(paste0(results_dir, 'day0_DE_genes_by_w5COCL2_presence.csv'))
de.cis <- read.csv(paste0(results_dir, 'day0_DE_genes_by_w5CIS_presence.csv'))

de.dabtram$labeling <- as.character(de.dabtram$labeling)
de.cocl2$labeling <- as.character(de.cocl2$labeling)
de.cis$labeling <- as.character(de.cis$labeling)
# ==============================================================================
# Plot
# ==============================================================================

p.dabtram <- ggplot(de.dabtram, aes(x = difference, y = log10pval)) + 
  geom_point() +
  geom_point(data = subset(de.dabtram, labeling %in% c("2", "3")), 
                               aes(color = labeling), size = 2) +
  scale_color_manual(values = c("1" = "black", "2" = "red", "3" = "blue")) +
  ggrepel::geom_text_repel(data = subset(de.dabtram, labeling %in% c("2", "3")),
                                    aes(label = name, color = labeling),
                                    box.padding = unit(0.5, 'lines'),
                                    point.padding = unit(1.6, 'lines'),
                                    size = 4,
                                    max.overlaps = 15) + 
  geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
                               color = "gray", linewidth=2) +
  geom_vline(xintercept=0, linetype="dashed", 
                               color = "gray", linewidth=2) +
  xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + 
  ylab("p-value (-Log10)") +
  theme_bw() + 
  Seurat::NoLegend()
ggsave(paste0(figure_dir, 'Supp_day0_DE_genes_by_w5DABTRAM_presence.png'), p.dabtram,
       width = 3, height = 3, dpi = 600)

p.cocl2 <- ggplot(de.cocl2, aes(x = difference, y = log10pval)) + 
  geom_point() +
  geom_point(data = subset(de.cocl2, labeling %in% c("2", "3")), 
             aes(color = labeling), size = 2) +
  scale_color_manual(values = c("1" = "black", "2" = "red", "3" = "blue")) +
  ggrepel::geom_text_repel(data = subset(de.cocl2, labeling %in% c("2", "3")),
                           aes(label = name, color = labeling),
                           box.padding = unit(0.5, 'lines'),
                           point.padding = unit(1.6, 'lines'),
                           max.overlaps = 50) + 
  geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
             color = "gray", linewidth=2) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "gray", linewidth=2) +
  xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + 
  ylab("p-value (-Log10)") +
  theme_bw() + 
  Seurat::NoLegend()
ggsave(paste0(figure_dir, 'Supp_day0_DE_genes_by_w5COCL2_presence.png'), p.cocl2,
       width = 3, height = 3, dpi = 600)


p.cis <- ggplot(de.cis, aes(x = difference, y = log10pval)) + 
  geom_point() +
  geom_point(data = subset(de.cis, labeling %in% c("2", "3")), 
             aes(color = labeling), size = 2) +
  scale_color_manual(values = c("1" = "black", "2" = "red", "3" = "blue")) +
  ggrepel::geom_text_repel(data = subset(de.cis, labeling %in% c("2", "3")),
                           aes(label = name, color = labeling),
                           box.padding = unit(0.5, 'lines'),
                           point.padding = unit(1.6, 'lines'),
                           max.overlaps = 50) + 
  geom_hline(yintercept=min(logpval_vec[which(pvalue_vec <= 0.05)]), linetype="dashed", 
             color = "gray", linewidth=2) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "gray", linewidth=2) +
  xlab("Mean difference in SAVER gene exp.: (Best-Worst)") + 
  ylab("p-value (-Log10)") +
  theme_bw() + 
  Seurat::NoLegend()
ggsave(paste0(figure_dir, 'Supp_day0_DE_genes_by_w5CIS_presence.png'), p.cis,
       width = 3, height = 3, dpi = 600)

