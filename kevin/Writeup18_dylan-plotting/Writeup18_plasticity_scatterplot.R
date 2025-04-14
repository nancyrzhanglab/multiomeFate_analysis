# https://github.com/nancyrzhanglab/multiomeFate_analysis/blob/emilia/emilia/Final/Fig2/Plot_lin_var_scores.R#L199
rm(list = ls())

set.seed(123)

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(grid)
library(ggthemes)
library(Seurat)
library(multiomeFate)

results_dir <- '~/project/Multiome_fate/out/emilia/task0_explore_lineage_variability_V2/'
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"


theme_Publication<- function(base_size=12, base_family="sans") {
  
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           plot.subtitle = element_text(face = "bold", hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold"),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(),
           axis.line.x = element_line(colour="black"),
           axis.line.y = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="white"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
           legend.direction = "horizontal",
           legend.key.size= unit(0.5, "cm"),
           legend.spacing = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
           strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
           strip.text = element_text(face="bold")
   ))
}


# read lineage variability scores
lin_var.day10_COCL2 <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_day10_COCL2_saver_sample.csv'))
lin_var.day10_COCL2_Shuffled <- read.csv(paste0(results_dir, 'day10_COCL2/lineage_variability_shuffledday10_COCL2_saver_sample.csv'))
lin_var.day10_COCL2_atac <- read.csv(paste0(results_dir, 'day10_COCL2/embeding_treatment/lineage_variability_day10_COCL2_peakvi.csv'))

# get lineage variability by RNA
lin_var.rna <- lin_var.day10_COCL2[, c('assigned_lineage', 'normalized_avg_eud_dist_by_shuffle')]
lin_var.rna$modality <- 'RNA'
colnames(lin_var.rna) <- c('assigned_lineage', 'lin_var.RNA', 'modality')

# get lineage variability by ATAC
lin_var.atac <- lin_var.day10_COCL2_atac[, c('assigned_lineage', 'normalized_avg_eud_dist_by_shuffle', 'n_cells')]
lin_var.atac$modality <- 'ATAC'
colnames(lin_var.atac) <- c('assigned_lineage', 'lin_var.ATAC', 'n_cells', 'modality')

# merge
rna.atac.comp <- merge(lin_var.rna, lin_var.atac, by = 'assigned_lineage')

################

p1 <- ggplot(rna.atac.comp, aes(x = lin_var.RNA, y = lin_var.ATAC)) +
  geom_point(aes(size = n_cells), shape = 21, fill = '#6DC49C') +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  stat_cor() +
  xlab('Lineage variability (RNA)') +
  ylab('Lineage variability (ATAC)') +
  theme_Publication()
ggsave(p1, file = paste0(plot_folder, 'Writuep18_plasticity_scatterplot.pdf'), width = 3.5, height = 4)


p1 <- ggplot(rna.atac.comp, aes(x = lin_var.RNA, y = lin_var.ATAC)) +
  geom_point(shape = 21, fill = '#6DC49C') +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  xlab('Lineage variability (RNA)') +
  ylab('Lineage variability (ATAC)') +
  theme_classic(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "gray80", linetype = "dotted", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(1,1,1,1)
  ) 
ggsave(p1, file = paste0(plot_folder, 'Writuep18_plasticity_scatterplot_cleaned.png'), width = 2, height = 2)

######################

# overlay with UMAP
all_data <- multiomeFate:::data_loader(which_files = c("rna", "fasttopics"))
umap <- as.data.frame(all_data[["ft.COCL2.umap"]]@cell.embeddings)

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)
metadat.day10_COCL2 <- metadat[metadat$dataset == 'day10_COCL2', ]
metadat.COCL2 <- metadat[metadat$dataset %in% c('day0', 'day10_COCL2', 'week5_COCL2'), ]

umap.COCL2 <- umap[metadat.COCL2$cell_id, ]
umap.day10_COCL2 <- umap[metadat.day10_COCL2$cell_id, ]

lin1.cells <- metadat.COCL2[metadat.COCL2$assigned_lineage == 'Lin70618', 'cell_id']
umap.day10_COCL2.lin1 <- umap.day10_COCL2[lin1.cells, ] # high variance

lin2.cells <- metadat.COCL2[metadat.COCL2$assigned_lineage == 'Lin34625', 'cell_id']
umap.day10_COCL2.lin2 <- umap.day10_COCL2[lin2.cells, ] # low variance

p2 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = 'gray90') +
  geom_point(data = umap.day10_COCL2.lin1,  size = 2, color = '#6DC49C') +
  ggtitle('High var.') +
  theme_Publication()

p3 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = 'gray90') +
  geom_point(data = umap.day10_COCL2.lin2,  size = 2, color = '#6DC49C') +
  ggtitle('Low var.') +
  theme_Publication()

ggsave(p2, file = paste0(plot_folder, 'Writuep18_plasticity_high-variablity_umap.pdf'), width = 3.5, height = 4)
ggsave(p3, file = paste0(plot_folder, 'Writuep18_plasticity_low-variablity_umap.pdf'), width = 3.5, height = 4)

p2 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = 'gray90') +
  geom_point(data = umap.day10_COCL2.lin1,  size = 2, color = '#6DC49C') +
  ggtitle('') + theme_classic(base_size = 8) 
p2 <- p2 + Seurat::NoAxes() + Seurat::NoLegend() + ggplot2::ggtitle("")
ggsave(p2, file = paste0(plot_folder, 'Writuep18_plasticity_high-variablity_umap_cleaned.png'), width = 2, height = 2)

p2 <- ggplot(umap.COCL2, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = 'gray90') +
  geom_point(data = umap.day10_COCL2.lin2,  size = 2, color = '#6DC49C') +
  ggtitle('') + theme_classic(base_size = 8) 
p2 <- p2 + Seurat::NoAxes() + Seurat::NoLegend() + ggplot2::ggtitle("")
ggsave(p2, file = paste0(plot_folder, 'Writuep18_plasticity_low-variablity_umap_cleaned.png'), width = 2, height = 2)

