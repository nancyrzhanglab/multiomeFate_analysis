rm(list = ls())

library(tidyverse)
library(ggplot2)


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
            # axis.ticks = element_blank(),
            # axis.line.x = element_line(colour="black"),
            # axis.line.y = element_line(colour="black"),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

use.pal_SCP_Signature1 <- c("#6BAED6", "#08519C", "#74C476", "#FCBBA1", "#EF3B2C", "#67000D", "#A63603", "grey75", "grey50")
names(use.pal_SCP_Signature1) <- c("SCP-Sensitive", "SCP-Sensitive2", "SCP 11", "SCP-Resistant1", "SCP-Resistant2", "SCP 11/SCP-Resistant2", "SCP-Resistant2;IFNmem-hi", "Plastic", "Missing")

use.pal_SCP_Signature2 <- c("#6BAED6", "#74C476", "#FCBBA1", "#EF3B2C", "#67000D", "grey75", "grey50")
names(use.pal_SCP_Signature2) <- c("SCP-Sen2", "SCP 11", "SCP-Res1", "SCP-Res2", "SCP 11/SCP-Res2", "Plastic", "Missing")

# ==============================================================================
# read data
# ==============================================================================
data_dir <- '/Users/emiliac/Dropbox/Epi Evolution Paper/Data/Lineage Trace/'
out_dir <- '~/Downloads/'

cell_anno <- readRDS(paste0(data_dir, 'Cell Barcode Metadata with Lineage Barcode Annotations 061124.rds'))
metadat.LB_ALL_FILT_Annotated <- readRDS(paste0(data_dir, 'metadat.LB_ALL_FILT_Annotated.rds'))

fp.nt <- readRDS(paste0(out_dir, 'final_fit_d14_d19_w_d14_size.rds'))
fp.jaki <- readRDS(paste0(out_dir, 'final_fit_d14_d19_JAKi_TBK1i_w_d14_size_scaled.fatepotential.rds'))

# ==============================================================================
# wrangle
# ==============================================================================

cell_anno <- cell_anno[cell_anno$In_Vivo == 'In Vivo', ]
cell_anno <- merge(cell_anno[, c('ID_Corrected', 'SEACell_Annotation_V3')], metadat.LB_ALL_FILT_Annotated, by = 'ID_Corrected')
cell_anno[, c('Label', 'LB_ID_short')] <- str_split_fixed(cell_anno$LB_ID, '_', 2)
cell_anno.d14.nt <- cell_anno[cell_anno$Label == 'Day 14 NT', ]
cell_anno.d14.jaki <- cell_anno[cell_anno$Label == 'Day 14 JAKi TBK1i', ]

# NT
fp.nt <- as.data.frame(fp.nt)
colnames(fp.nt) <- 'FatePotential.d14.d19.NT'
fp.nt$ID_Corrected <- rownames(fp.nt)

fp.nt <- merge(fp.nt, cell_anno, by = 'ID_Corrected')

# JAKi TBK1i
fp.jaki <- as.data.frame(fp.jaki)
colnames(fp.jaki) <- 'FatePotential.d14.d19.JAKi'
fp.jaki$ID_Corrected <- rownames(fp.jaki)

fp.jaki <- merge(fp.jaki, cell_anno, by = 'ID_Corrected')

# future lineage size
lin.size <- metadat.LB_ALL_FILT_Annotated %>% 
  group_by(LB_ID, Label) %>%
  summarise(Size = n())

lin.size.d14.nt <- lin.size[lin.size$Label == 'Day 14 NT', ]
lin.size.d14.nt <- lin.size.d14.nt[order(lin.size.d14.nt$Size, decreasing = TRUE), ]
lin.size.d14.jaki <- lin.size[lin.size$Label == 'Day 14 JAKi TBK1i', ]
lin.size.d14.jaki <- lin.size.d14.jaki[order(lin.size.d14.jaki$Size, decreasing = TRUE), ]

lin.size.d19.nt <- lin.size[lin.size$Label == 'Day 19 NT', ]
lin.size.d19.jaki <- lin.size[lin.size$Label == 'Day 19 JAKi TBK1i', ]

colnames(lin.size.d19.nt) <- c('LB_ID', 'Label', 'Lin.Size.d19.NT')
colnames(lin.size.d19.jaki) <- c('LB_ID', 'Label', 'Lin.Size.d19.JAKi')

lin.size.d19.nt[, c('Label', 'LB_ID_short')] <- str_split_fixed(lin.size.d19.nt$LB_ID, '_', 2)
lin.size.d19.jaki[, c('Label', 'LB_ID_short')] <- str_split_fixed(lin.size.d19.jaki$LB_ID, '_', 2)

lin.size.d19.nt <- merge(lin.size.d19.nt[, c('LB_ID_short', 'Lin.Size.d19.NT')], cell_anno.d14.nt, by = c('LB_ID_short'))
lin.size.d19.jaki <- merge(lin.size.d19.jaki[, c('LB_ID_short', 'Lin.Size.d19.JAKi')], cell_anno.d14.jaki, by = c('LB_ID_short'))

# ==============================================================================
# plot 
# ==============================================================================

ggplot(metadat.LB_ALL_FILT_Annotated, aes(x = UMAP1_postHarmony, y = UMAP2_postHarmony, color = Day)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#FF9D23", "#493D9E")) +
  xlim(-5, 7) +
  ylim(-3, 8) +
  theme_Publication()

ggplot(metadat.LB_ALL_FILT_Annotated, aes(x = UMAP1_postHarmony, y = UMAP2_postHarmony, color = Day)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = c("#FF9D23", "#493D9E")) +
  xlim(-5, 7) +
  ylim(-3, 8) +
  theme_Publication()


ggplot(cell_anno, aes(x = UMAP1_postHarmony, y = UMAP2_postHarmony, color = SEACell_Annotation_V3)) + 
  geom_point(size = 0.5) +
  scale_color_manual(values = use.pal_SCP_Signature2) +
  xlim(-5, 7) +
  ylim(-3, 8) +
  theme_Publication()

fp.nt <- fp.nt[order(fp.nt$FatePotential.d14.d19.NT), ]
ggplot(fp.nt, aes(x = UMAP1_postHarmony, y = UMAP2_postHarmony, color = FatePotential.d14.d19.NT)) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "bisque", high = "red", midpoint = -1.6) +
  xlim(-5, 7) +
  ylim(-3, 8) +
  theme_Publication()

fp.jaki <- fp.jaki[order(fp.jaki$FatePotential.d14.d19.JAKi), ]
ggplot(fp.jaki, aes(x = UMAP1_postHarmony, y = UMAP2_postHarmony, color = FatePotential.d14.d19.JAKi)) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "bisque", high = "red", midpoint = -1.4) +
  xlim(-5, 7) +
  ylim(-3, 8) +
  theme_Publication()


lin.size.large.d14.nt <- head(lin.size.d14.nt, 15)
fp.var.nt <- fp.nt %>% 
  group_by(LB_ID) %>%
  summarise(FP.var = var(FatePotential.d14.d19.NT),
            FP.median = median(FatePotential.d14.d19.NT)) %>% 
  arrange(desc(FP.var))
fp.var.large.d14.nt <- head(fp.var.nt, 15)
fp.var.large.d14.nt <- fp.var.large.d14.nt[order(fp.var.large.d14.nt$FP.median, decreasing = T), ]

fp.nt$lin.plot <- ifelse(fp.nt$LB_ID %in% fp.var.large.d14.nt$LB_ID, fp.nt$LB_ID, 'other')
fp.nt$lin.plot <- factor(fp.nt$lin.plot, levels = c(fp.var.large.d14.nt$LB_ID, 'other'))
ggplot(fp.nt, aes(x = lin.plot, y = FatePotential.d14.d19.NT)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.5, outlier.shape = NA, size = 1) +
  geom_jitter(width = 0.2, size = 0.5, color = 'gray') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5))
# ggsave('~/Downloads/FP.var.large.d14.nt.png', width = 8, height = 7, dpi = 300)


fp.var.jaki <- fp.jaki %>% 
  group_by(LB_ID) %>%
  summarise(FP.var = var(FatePotential.d14.d19.JAKi),
            FP.median = median(FatePotential.d14.d19.JAKi)) %>% 
  arrange(desc(FP.var))
fp.var.large.d14.jaki <- head(fp.var.jaki, 15)
fp.var.large.d14.jaki <- fp.var.large.d14.jaki[order(fp.var.large.d14.jaki$FP.median, decreasing = T), ]

fp.jaki$lin.plot <- ifelse(fp.jaki$LB_ID %in% fp.var.large.d14.jaki$LB_ID, fp.jaki$LB_ID, 'other')
fp.jaki$lin.plot <- factor(fp.jaki$lin.plot, levels = c(fp.var.large.d14.jaki$LB_ID, 'other'))
ggplot(fp.jaki, aes(x = lin.plot, y = FatePotential.d14.d19.JAKi)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.5, outlier.shape = NA, size = 1) +
  geom_jitter(width = 0.2, size = 0.5, color = 'gray') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5))
ggsave('~/Downloads/FP.var.large.d14.jaki.png', width = 8, height = 7.5, dpi = 300)

lin.size.d19.nt <- lin.size.d19.nt[order(lin.size.d19.nt$Lin.Size.d19.NT), ]
hist(lin.size.d19.nt$Lin.Size.d19, breaks = 50)
lin.size.d19.nt$Lin.Size.d19.NT.log10 <- log10(lin.size.d19.nt$Lin.Size.d19 + 1)
ggplot(lin.size.d19.nt, aes(x = UMAP1_postHarmony, y = UMAP2_postHarmony, color = Lin.Size.d19.NT.log10)) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "bisque", high = "red", midpoint = 1.041393) +
  xlim(-5, 7) +
  ylim(-3, 8) +
  theme_Publication()

lin.size.d19.jaki <- lin.size.d19.jaki[order(lin.size.d19.jaki$Lin.Size.d19.JAKi), ]
lin.size.d19.jaki$Lin.Size.d19.JAKi.log10 <- log10(lin.size.d19.jaki$Lin.Size.d19.JAKi + 1)
ggplot(lin.size.d19.jaki, aes(x = UMAP1_postHarmony, y = UMAP2_postHarmony, color = Lin.Size.d19.JAKi.log10)) +
  geom_point(size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "bisque", high = "red", midpoint = 0.698970) +
  xlim(-5, 7) +
  ylim(-3, 8) +
  theme_Publication()
