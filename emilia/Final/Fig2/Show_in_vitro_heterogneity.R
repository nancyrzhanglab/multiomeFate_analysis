rm(list = ls())

library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
# results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/Signatures/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig2/'

data_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/kevin/Writeup10a/'
results_dir <- '/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/emilia/task2_correlate_fate_potential_and_features_V2/'

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


remove_unassigned_cells <- TRUE

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir,'Writeup10a_data_wnn.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver
all_data[["wnn.umap"]] <- all_data_wnn

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

mitf_axl <- read.csv(paste0(ref_dir, 'Triosh_2016_MITF_AXL_Programs.csv'))

# ==============================================================================
# Wrangle data
# ==============================================================================
mitf_program <- unique(na.omit(mitf_axl$MITF_Program[mitf_axl$MITF_Program %in% all_data_saver@data@Dimnames[[1]]]))
axl_program <- unique(na.omit(mitf_axl$AXL_Program[mitf_axl$AXL_Program %in% all_data_saver@data@Dimnames[[1]]]))

# ==============================================================================
# Calcualte module scores 
# ==============================================================================
metadat <- all_data@meta.data
metadat$cell_barcode <- rownames(metadat)

scores <- ScoreSignatures_UCell(all_data@assays[["Saver"]]@data, 
                                features=c(list(mitf_program), 
                                           list(axl_program)))

colnames(scores) <- c('MITF_Program', 'AXL_Program')

# write.csv(scores, paste0(results_dir, 'MITF_AXL_UCell_scores.csv'))
# ==============================================================================
# Plot UCell scores
# ==============================================================================

scores <- read.csv(paste0(results_dir, 'MITF_AXL_UCell_scores.csv'), row.names = 1)
scores.df <- as.data.frame(scores)
scores.df$cell_barcode <- rownames(scores.df)

scores.df.w.metadat <- merge(scores.df, metadat[, c('cell_barcode', 'dataset', 'assigned_lineage',
                                                    'fatepotential_DABTRAM_d0_d10', 'fatepotential_COCL2_d0_d10', 'fatepotential_CIS_d0_d10',
                                                    'fatepotential_DABTRAM_d10_w5', 'fatepotential_COCL2_d10_w5', 'fatepotential_CIS_d10_w5')], by = 'cell_barcode')
scores.df.long <- melt(scores.df.w.metadat, id.vars = c('cell_barcode', 'dataset', 'assigned_lineage',
                                                        'fatepotential_DABTRAM_d0_d10', 'fatepotential_COCL2_d0_d10', 'fatepotential_CIS_d0_d10',
                                                        'fatepotential_DABTRAM_d10_w5', 'fatepotential_COCL2_d10_w5', 'fatepotential_CIS_d10_w5'))

umap <- as.data.frame(all_data@reductions[["wnn.umap"]]@cell.embeddings)
umap$cell_barcode <- rownames(umap)

scores.df.long <- merge(scores.df.long, umap, by = 'cell_barcode')

scores.df.long.AXL<- scores.df.long[scores.df.long$variable == 'AXL_Program', ]
scores.df.long.MITF<- scores.df.long[scores.df.long$variable == 'MITF_Program', ]

# Day10 dabtram AXL
scores.df.long.day10Dabtram <- scores.df.long.AXL[scores.df.long.AXL$dataset == 'day10_DABTRAM', ]
scores.df.long.day10Dabtram$value <- scale(scores.df.long.day10Dabtram$value)
scores.df.long.day10Dabtram.lin <- scores.df.long.day10Dabtram %>% 
  group_by(assigned_lineage) %>% 
  summarise(mean = mean(value),
            var = var(value),
            lin_size = n()) %>% 
  filter(lin_size > 10)
scores.df.long.day10Dabtram.lin <- scores.df.long.day10Dabtram.lin[order(scores.df.long.day10Dabtram.lin$var), ]
lins_to_plot <- c(head(scores.df.long.day10Dabtram.lin, 4)$assigned_lineage,
                  tail(scores.df.long.day10Dabtram.lin, 4)$assigned_lineage)

p1 <- ggplot(scores.df.long.day10Dabtram, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_Publication() +
  theme(panel.grid = element_blank(),
        legend.position = 'none')

scores.df.long.day10Dabtram.to_plot <- scores.df.long.day10Dabtram[scores.df.long.day10Dabtram$assigned_lineage %in% lins_to_plot, ]
p1.1 <- ggplot(scores.df.long.day10Dabtram.to_plot, aes(x = reorder(assigned_lineage, -value, median), y = value)) +
  geom_violin(scale = 'width')+
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  ylab('AXL Program') +
  theme_classic() +
  xlab('') +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  theme(strip.text = element_text(size = 6))


# Day10 COCL2 AXL
scores.df.long.day10COCL2 <- scores.df.long.AXL[scores.df.long.AXL$dataset == 'day10_COCL2', ]
scores.df.long.day10COCL2$value <- scale(scores.df.long.day10COCL2$value)

p2 <- ggplot(scores.df.long.day10COCL2, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_Publication() +
  theme(panel.grid = element_blank(),
        legend.position = 'none')

scores.df.long.day10COCL2.lin <- scores.df.long.day10COCL2 %>% 
  group_by(assigned_lineage) %>% 
  summarise(mean = mean(value),
            var = var(value),
            lin_size = n()) %>% 
  filter(lin_size > 10) %>% 
  drop_na()
scores.df.long.day10COCL2.lin <- scores.df.long.day10COCL2.lin[order(scores.df.long.day10COCL2.lin$var), ]
lins_to_plot <- c(head(scores.df.long.day10COCL2.lin, 4)$assigned_lineage,
                  tail(scores.df.long.day10COCL2.lin, 4)$assigned_lineage)

scores.df.long.day10COCL2.to_plot <- scores.df.long.day10COCL2[scores.df.long.day10COCL2$assigned_lineage %in% lins_to_plot, ]
p2.1 <- ggplot(scores.df.long.day10COCL2.to_plot, aes(x = reorder(assigned_lineage, -value, median), y = value)) +
  geom_violin(scale = 'width')+
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  ylab('AXL Program') +
  xlab('') +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  theme(strip.text = element_text(size = 6))


# Day10 CIS AXL
scores.df.long.day10CIS <- scores.df.long.AXL[scores.df.long.AXL$dataset == 'day10_CIS', ]
scores.df.long.day10CIS$value <- scale(scores.df.long.day10CIS$value)

p3 <- ggplot(scores.df.long.day10CIS, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_Publication() +
  theme(panel.grid = element_blank(),
        legend.position = 'none')

scores.df.long.day10CIS.lin <- scores.df.long.day10CIS %>% 
  group_by(assigned_lineage) %>% 
  summarise(mean = mean(value),
            var = var(value),
            lin_size = n()) %>% 
  filter(lin_size > 10)
scores.df.long.day10CIS.lin <- scores.df.long.day10CIS.lin[order(scores.df.long.day10CIS.lin$var), ]
lins_to_plot <- c(head(scores.df.long.day10CIS.lin, 4)$assigned_lineage,
                  tail(scores.df.long.day10CIS.lin, 4)$assigned_lineage)

scores.df.long.day10CIS.to_plot <- scores.df.long.day10CIS[scores.df.long.day10CIS$assigned_lineage %in% lins_to_plot, ]
p3.1 <- ggplot(scores.df.long.day10CIS.to_plot, aes(x = reorder(assigned_lineage, -value, median), y = value)) +
  geom_violin(scale = 'width')+
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  ylab('AXL Program') +
  xlab('') +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  theme(strip.text = element_text(size = 6))



p4 <- grid.arrange(p1, p2, p3, ncol = 1)
ggsave(paste0(figure_dir, 'MITF_AXL_D10.pdf'), p7, width = 2, height = 5)

p5 <- grid.arrange(p1.1, p2.1, p3.1, ncol = 1)
ggsave(paste0(figure_dir, 'MITF_AXL_D10_violin.pdf'), p5, width = 3, height = 5.5)

p6 <- ggarrange(p4, p5, ncol = 2, widths = c(0.5, 1))

ggsave(paste0(figure_dir, 'MITF_AXL_D10_umap_plus_violin.pdf'), p6, width = 6, height = 5.5)


