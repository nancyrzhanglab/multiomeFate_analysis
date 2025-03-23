rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggdensity)
library(ggplot2)
library(RColorBrewer)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'


remove_unassigned_cells <- TRUE
treatment <- 'DABTRAM'

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
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
# ==============================================================================
# Read data
# ==============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
metadat <- all_data@meta.data

# fate bias from d0 to week5 adapting
df.bias <- read.csv(paste0(out_dir, 'adapting_bias_thres_0_', treatment, '.csv'))

# gavish scores
rna.gavish.mp <- read.csv(paste0(results_dir, 'GavishMP_UCell_scores.csv'), row.names = 1)

# ==============================================================================
# Wrangle
# ==============================================================================
df.bias.high <- df.bias %>% filter(bias > 0.5)
df.bias.low <- df.bias %>% filter(bias < 0.5)

rna.gavish.mp.day0.bias.high <- rna.gavish.mp[df.bias.high$cell_id, ] 
rna.gavish.mp.day0.bias.low <- rna.gavish.mp[df.bias.low$cell_id, ]

metadat$cell_id <- rownames(metadat)
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]
metadat.week5_DABTRAM <- metadat[metadat$dataset == 'week5_DABTRAM', ]
rna.gavish.mp.day10_DABTRAM <- rna.gavish.mp[metadat.day10_DABTRAM$cell_id, ]
rna.gavish.mp.week5_DABTRAM <- rna.gavish.mp[metadat.week5_DABTRAM$cell_id, ]

rna.gavish.mp.day0.bias.high$category <- 'day0.bias.high'
rna.gavish.mp.day0.bias.low$category <- 'day0.bias.low'
rna.gavish.mp.day10_DABTRAM$category <- 'day10_DABTRAM'
rna.gavish.mp.week5_DABTRAM$category <- 'week5_DABTRAM'

to_plot <- rbind(rna.gavish.mp.day0.bias.high, rna.gavish.mp.day0.bias.low, 
                 rna.gavish.mp.day10_DABTRAM, rna.gavish.mp.week5_DABTRAM)
to_plot$category <- factor(to_plot$category, levels = c('day0.bias.low', 'day0.bias.high', 'day10_DABTRAM', 'week5_DABTRAM'))
# ==============================================================================
# Plot
# ==============================================================================
p1 <- ggplot(to_plot, aes(x = category, y = `Interferon.MHC.II..I._UCell`)) +
        geom_violin(scale = 'width') +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        theme_Publication()

p2 <- ggplot(to_plot, aes(x = category, y = `Skin.pigmentation_UCell`)) +
        geom_violin(scale = 'width') +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        theme_Publication()

ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = TRUE, legend = 'bottom')
