library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '~/Downloads/'

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

gavish.mp <- read_csv(paste0(ref_dir, 'Gavish_Malignant_Meta_Programs.csv'))
isg_rs <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Gavish Clinical Analysis/Resources/ISG.RS.txt', col_names = F)
isg_mem <- read_table('/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Persistent IFN ISG Groups/Memory ISGs Human.csv', col_names = F)

colnames(isg_rs) <- 'Gene'
colnames(isg_mem) <- 'Gene'

mitf_axl <- read.csv('~/Downloads/Triosh_2016_MITF_AXL_Programs.csv')
# ==============================================================================
# Wrangle data
# ==============================================================================
gavish.mp.list <- list()
for (mp in colnames(gavish.mp)) {
  gavish.mp.list[[mp]] <- c(gavish.mp[[mp]] %>% as.character())
}

gavish.mp.list <- lapply(1:length(gavish.mp.list), function(x) {
  gs <- gavish.mp.list[[x]]
  gs <- unique(na.omit(gs[gs %in% all_data_saver@data@Dimnames[[1]]]))
  return(gs)
})
names(gavish.mp.list) <- colnames(gavish.mp)

isg_rs <- unique(na.omit(isg_rs$Gene[isg_rs$Gene %in% all_data_saver@data@Dimnames[[1]]]))
isg_mem <- unique(na.omit(isg_mem$Gene[isg_mem$Gene %in% all_data_saver@data@Dimnames[[1]]]))

mitf_program <- unique(na.omit(mitf_axl$MITF_Program[mitf_axl$MITF_Program %in% all_data_saver@data@Dimnames[[1]]]))
axl_program <- unique(na.omit(mitf_axl$AXL_Program[mitf_axl$AXL_Program %in% all_data_saver@data@Dimnames[[1]]]))
# ==============================================================================
# Calcualte module scores 
# ==============================================================================
all_data1 <- subset(all_data, dataset == 'day0')
all_data1 <- AddModuleScore(all_data1, 
                           features = gavish.mp.list, 
                           assay = "Saver",
                           name = names(gavish.mp.list))

all_data1 <- AddModuleScore(all_data1, 
                            features = list(isg_rs), 
                            assay = "Saver",
                            name = 'ISG.RS')

all_data1 <- AddModuleScore(all_data1, 
                            features = list(isg_mem), 
                            assay = "Saver",
                            name = 'ISG.Mem')

all_data1 <- AddModuleScore(all_data1, 
                            features = list(mitf_program), 
                            assay = "Saver",
                            name = 'MITF_program')

all_data1 <- AddModuleScore(all_data1, 
                            features = list(axl_program), 
                            assay = "Saver",
                            name = 'AXL_program')

metadat <- all_data@meta.data
scores <- ScoreSignatures_UCell(all_data@assays[["Saver"]]@data, 
                                features=c(gavish.mp.list, 
                                           list(isg_rs), 
                                           list(isg_mem), 
                                           list(mitf_program), 
                                           list(axl_program)))

colnames(scores) <- c(names(gavish.mp.list), 'ISG.RS', 'ISG.Mem', 'MITF_Program', 'AXL_Program')

# ==============================================================================
# Plot module scores
# ==============================================================================
metadat <- metadat[, c('dataset', paste0(colnames(gavish.mp), seq(1, 41)), 'ISG.RS1', 'ISG.Mem1', 'MITF_program1', 'AXL_program1',
                       'fatepotential_DABTRAM_d0_d10', 'fatepotential_COCL2_d0_d10', 'fatepotential_CIS_d0_d10')]
metadat$cell_barcode <- rownames(metadat)
metadat_melt <- melt(metadat, id.vars = c('cell_barcode', 'dataset', 'fatepotential_DABTRAM_d0_d10', 'fatepotential_COCL2_d0_d10', 'fatepotential_CIS_d0_d10'))

ggplot(metadat_melt, aes(x = dataset, y = value)) +
  geom_violin(scale = 'width') +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  facet_wrap(. ~ variable, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

val_to_plot <- c('Cell Cycle - G2/M1', 'Cell Cycle - G1/S2', 
                 'Cell Cycle HMG-rich3', 'Stress5', 'Hypoxia6', 'Stress (in vitro)7',
                 'EMT-I12', 'EMT-II13', 'EMT-III14', 'EMT-IV15', 
                 'Interferon/MHC-II (I)17', 'Interferon/MHC-II (II)18', 'Epithelial Senescence19')
metadat_melt <- metadat_melt[metadat_melt$variable %in% val_to_plot, ]
ggplot(metadat_melt, aes(x = variable, y = value)) +
  geom_violin(scale = 'width') +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(metadat_melt, aes(x = value, y = fatepotential_COCL2_d0_d10)) +
  geom_point(size = 1) +
  stat_cor(method = 'spearman', label.x.npc = "left", label.y.npc = "bottom") +
  geom_smooth(method = 'lm') +
  facet_wrap(. ~ variable, scales = 'free_x') +
  theme_bw()

FeaturePlot(all_data1, 
            features = c(val_to_plot), 
            order = T,
            reduction = "wnn.umap")

FeaturePlot(all_data1, 
            features = c('AXL_program1', 'MITF_program1', 'fatepotential_DABTRAM_d0_d10'), 
            reduction = "wnn.umap")

# ==============================================================================
# Plot UCell scores
# ==============================================================================
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
scores.df.long <- scores.df.long[scores.df.long$variable %in% c('MITF_Program', 'AXL_Program', 
                                                                'Stress', 'Hyposia', 'EMT-I', 'EMT-II', 'EMT-III', 'EMT-IV', 'Interferon/MHC-II (I)', 'Interferon/MHC-II (II)',
                                                                'ISG.RS'), ]
ggplot(scores.df.long, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  facet_wrap(. ~ variable, scales = 'free') +
  theme_bw() +
  theme(panel.grid = element_blank())

scores.df.long.AXL<- scores.df.long[scores.df.long$variable == 'AXL_Program', ]
scores.df.long.MITF<- scores.df.long[scores.df.long$variable == 'MITF_Program', ]

# scores.df.long.AXL$value_scale <- scale(scores.df.long.AXL$value)
# ggplot(scores.df.long.AXL, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
#   geom_point(aes(color = value_scale), size = 1) +
#   scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
#                         limits = c(-2, 2),
#                         oob = scales::squish) +
#   facet_wrap(. ~ dataset) +
#   theme_classic() +
#   theme(panel.grid = element_blank(),
#         legend.position = 'none')

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
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = 'none')

scores.df.long.day10Dabtram.to_plot <- scores.df.long.day10Dabtram[scores.df.long.day10Dabtram$assigned_lineage %in% lins_to_plot, ]
p1.1 <- ggplot(scores.df.long.day10Dabtram.to_plot, aes(x = reorder(assigned_lineage, -value, median), y = value)) +
  geom_violin(scale = 'width')+
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  ylab('AXL Program') +
  theme_classic() +
  xlab('') +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  theme(strip.text = element_text(size = 6))

# Day10 dabtram MITF
scores.df.long.day10Dabtram <- scores.df.long.MITF[scores.df.long.MITF$dataset == 'day10_DABTRAM', ]
scores.df.long.day10Dabtram$value <- scale(scores.df.long.day10Dabtram$value)

p2 <- ggplot(scores.df.long.day10Dabtram, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = 'none')

# Day10 COCL2 AXL
scores.df.long.day10COCL2 <- scores.df.long.AXL[scores.df.long.AXL$dataset == 'day10_COCL2', ]
scores.df.long.day10COCL2$value <- scale(scores.df.long.day10COCL2$value)

p3 <- ggplot(scores.df.long.day10COCL2, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_classic() +
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
  geom_jitter(width = 0.1, alpha = 0.2) +
  ylab('AXL Program') +
  xlab('') +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  theme(strip.text = element_text(size = 6))


# Day10 COCL2 MITF
scores.df.long.day10COCL2 <- scores.df.long.MITF[scores.df.long.MITF$dataset == 'day10_COCL2', ]
scores.df.long.day10COCL2$value <- scale(scores.df.long.day10COCL2$value)

p4 <- ggplot(scores.df.long.day10COCL2, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = 'none')


# Day10 CIS AXL
scores.df.long.day10CIS <- scores.df.long.AXL[scores.df.long.AXL$dataset == 'day10_CIS', ]
scores.df.long.day10CIS$value <- scale(scores.df.long.day10CIS$value)

p5 <- ggplot(scores.df.long.day10CIS, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_classic() +
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
  geom_jitter(width = 0.1, alpha = 0.2) +
  ylab('AXL Program') +
  xlab('') +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  theme(strip.text = element_text(size = 6))

# Day10 CIS MITF
scores.df.long.day10CIS <- scores.df.long.MITF[scores.df.long.MITF$dataset == 'day10_CIS', ]
scores.df.long.day10CIS$value <- scale(scores.df.long.day10CIS$value)

p6 <- ggplot(scores.df.long.day10CIS, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = value), size = 1) +
  scale_color_gradient2(low = '#604CC3', mid = 'gray', high = '#FFA500',
                        limits = c(-2, 2),
                        oob = scales::squish) +
  theme_classic() +
  theme(panel.grid = element_blank())


p7 <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
# ggsave('~/Downloads/MITF_AXL_D10.png', p7, width = 4, height = 5, dpi = 300)

p8 <- grid.arrange(p1.1, p2.1, p3.1, ncol = 1)
ggsave('~/Downloads/MITF_AXL_D10_violin.png', p8, width = 3, height = 5.5, dpi = 300)

ggplot(scores.df.long, aes(x = dataset, y = value)) +
  geom_violin(scale = 'width') +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  facet_wrap(. ~ variable, scales = 'free_y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(scores.df.long, aes(x = fatepotential_CIS_d10_w5, y = value)) +
  geom_point(size = 1) +
  geom_smooth(method = 'lm') +
  # stat_cor(method = 'spearman', label.x.npc = "left", label.y.npc = "bottom") +
  facet_wrap(. ~ variable, scales = 'free_y') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))

scores.df.w.metadat.d0 <- scores.df.w.metadat[scores.df.w.metadat$dataset == 'day0', ]
scores.df.w.metadat.d0 <- merge(scores.df.w.metadat.d0, umap, by = 'cell_barcode')
scores.df.w.metadat.d0$AXL_Program_scaled <- scale(scores.df.w.metadat.d0$AXL_Program)
scores.df.w.metadat.d0$MITF_Program_scaled <- scale(scores.df.w.metadat.d0$MITF_Program)
ggplot(scores.df.w.metadat.d0, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = fatepotential_COCL2_d0_d10), size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(scores.df.w.metadat.d0, aes(x = MITF_Program_scaled, y = AXL_Program_scaled)) +
  geom_point(aes(color = fatepotential_DABTRAM_d0_d10), size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank())



scores.df.w.metadat.d10_DABTRAM <- scores.df.w.metadat[scores.df.w.metadat$dataset == 'day10_DABTRAM', ]
ggplot(scores.df.w.metadat.d10_DABTRAM, aes(x = MITF_Program, y = AXL_Program)) +
  geom_point(aes(color = fatepotential_DABTRAM_d10_w5), size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(scores.df.w.metadat.d10_DABTRAM, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = ISG_RS_scaled), size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank())

scores.df.w.metadat.d10_DABTRAM <- merge(scores.df.w.metadat.d10_DABTRAM, umap, by = 'cell_barcode')
scores.df.w.metadat.d10_DABTRAM$AXL_Program_scaled <- scale(scores.df.w.metadat.d10_DABTRAM$AXL_Program)
scores.df.w.metadat.d10_DABTRAM$MITF_Program_scaled <- scale(scores.df.w.metadat.d10_DABTRAM$MITF_Program)
scores.df.w.metadat.d10_DABTRAM$ISG_RS_scaled <- scale(scores.df.w.metadat.d10_DABTRAM$ISG.RS)

ggplot(scores.df.w.metadat.d10_DABTRAM, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = MITF_Program_scaled), size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(scores.df.w.metadat.d10_DABTRAM, aes(x = fatepotential_DABTRAM_d10_w5, y = ISG_RS_scaled)) +
  geom_point(size = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())


ggplot(scores.df.w.metadat, aes(x = fatepotential_DABTRAM_d10_w5, y = variable)) +
  geom_point(size = 1) +
  stat_cor(method = 'spearman', label.x.npc = "left", label.y.npc = "bottom") +
  theme_bw() +
  theme(panel.grid = element_blank())

scores.df.w.metadat.DABTRAM <- scores.df.w.metadat[scores.df.w.metadat$dataset %in% c('day0', 'day10_DABTRAM', 'week5_DABTRAM'), ]
scores.df.w.metadat.DABTRAM <- merge(scores.df.w.metadat.DABTRAM, umap, by = 'cell_barcode')
scores.df.w.metadat.DABTRAM$ISG_RS_scaled <- scale(scores.df.w.metadat.DABTRAM$ISG.RS)
scores.df.w.metadat.DABTRAM$Stress_scaled <- scale(scores.df.w.metadat.DABTRAM$Stress)
scores.df.w.metadat.DABTRAM$AXL_Program_scaled <- scale(scores.df.w.metadat.DABTRAM$AXL_Program)
scores.df.w.metadat.DABTRAM$MITF_Program_scaled <- scale(scores.df.w.metadat.DABTRAM$MITF_Program)
scores.df.w.metadat.DABTRAM$`Interferon/MHC-II (I)_scaled` <- scale(scores.df.w.metadat.DABTRAM$`Interferon/MHC-II (I)`)
scores.df.w.metadat.DABTRAM$`Interferon/MHC-II (II)_scaled` <- scale(scores.df.w.metadat.DABTRAM$`Interferon/MHC-II (II)`)

ggplot(scores.df.w.metadat.DABTRAM, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = `Stress_scaled`), size = 1) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  facet_wrap(. ~ dataset) +
  theme_bw() +
  theme(panel.grid = element_blank())


scores.df.w.metadat.COCL2 <- scores.df.w.metadat[scores.df.w.metadat$dataset %in% c('day0', 'day10_COCL2', 'week5_COCL2'), ]
scores.df.w.metadat.COCL2 <- merge(scores.df.w.metadat.COCL2, umap, by = 'cell_barcode')
scores.df.w.metadat.COCL2$ISG_RS_scaled <- scale(scores.df.w.metadat.COCL2$ISG.RS)
scores.df.w.metadat.COCL2$Stress_scaled <- scale(scores.df.w.metadat.COCL2$Stress)
scores.df.w.metadat.COCL2$AXL_Program_scaled <- scale(scores.df.w.metadat.COCL2$AXL_Program)
scores.df.w.metadat.COCL2$MITF_Program_scaled <- scale(scores.df.w.metadat.COCL2$MITF_Program)
scores.df.w.metadat.COCL2$`Interferon/MHC-II (I)_scaled` <- scale(scores.df.w.metadat.COCL2$`Interferon/MHC-II (I)`)
scores.df.w.metadat.COCL2$`Interferon/MHC-II (II)_scaled` <- scale(scores.df.w.metadat.COCL2$`Interferon/MHC-II (II)`)
scores.df.w.metadat.COCL2$`EMT-I_scaled` <- scale(scores.df.w.metadat.COCL2$`EMT-I`)

ggplot(scores.df.w.metadat.COCL2, aes(x = wnnUMAP_1, y = wnnUMAP_2)) +
  geom_point(aes(color = `EMT-I_scaled`), size = 0.5) +
  scale_color_gradient2(low = 'blue', mid = 'gray', high = 'red') +
  facet_wrap(. ~ dataset) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot(scores.df.w.metadat.COCL2, aes(x = fatepotential_COCL2_d10_w5, y = `EMT-I_scaled`)) +
  geom_point(size = 1) +
  stat_cor() +
  theme_bw() +
  theme(panel.grid = element_blank())
