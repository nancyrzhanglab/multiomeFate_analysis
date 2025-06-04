rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
out_dir1 <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
out_dir2 <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

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

# =============================================================================
# Read and wrangle data
# =============================================================================

gavish.mp.scores <- read.csv(paste0(out_dir1, 'GavishMP_UCell_scores.csv'), row.names = 1)
df.bias.DABTRAM <- read.csv(paste0(out_dir2, 'adapting_bias_thres_0_DABTRAM.csv'))
df.bias.COCL2 <- read.csv(paste0(out_dir2, 'adapting_bias_thres_0_COCL2.csv'))
df.bias.CIS <- read.csv(paste0(out_dir2, 'adapting_bias_thres_0_CIS.csv'))

gavish.mps <- colnames(gavish.mp.scores)
gavish.mp.scores$cell_id <- rownames(gavish.mp.scores)
df.bias.DABTRAM <- merge(df.bias.DABTRAM, gavish.mp.scores, by = 'cell_id')
df.bias.COCL2 <- merge(df.bias.COCL2, gavish.mp.scores, by = 'cell_id')
df.bias.CIS <- merge(df.bias.CIS, gavish.mp.scores, by = 'cell_id')

# ==============================================================================
# Correlation test
# ==============================================================================

# DABTRAM
cor.vec.dabtram <- sapply(gavish.mps, function(x) {
  res <- cor.test(df.bias.DABTRAM[, x], df.bias.DABTRAM$bias, method = 'spearman')
  return(c(res$estimate, res$p.value))
})
cor.vec.dabtram <- as.data.frame(t(cor.vec.dabtram))
rownames(cor.vec.dabtram) <- gsub('_UCell', '.rho', gavish.mps)
colnames(cor.vec.dabtram) <- c('cor.DABTRAM', 'pval.DABTRAM')
cor.vec.dabtram$MetaProgram <- rownames(cor.vec.dabtram)

# COCL2
cor.vec.cocl2 <- sapply(gavish.mps, function(x) {
  res <- cor.test(df.bias.COCL2[, x], df.bias.COCL2$bias, method = 'spearman')
  return(c(res$estimate, res$p.value))
})
cor.vec.cocl2 <- as.data.frame(t(cor.vec.cocl2))
rownames(cor.vec.cocl2) <- gsub('_UCell', '.rho', gavish.mps)
colnames(cor.vec.cocl2) <- c('cor.COCL2', 'pval.COCL2')
cor.vec.cocl2$MetaProgram <- rownames(cor.vec.cocl2)

# CIS
cor.vec.cis <- sapply(gavish.mps, function(x) {
  res <- cor.test(df.bias.CIS[, x], df.bias.CIS$bias, method = 'spearman')
  return(c(res$estimate, res$p.value))
})
cor.vec.cis <- as.data.frame(t(cor.vec.cis))
rownames(cor.vec.cis) <- gsub('_UCell', '.rho', gavish.mps)
colnames(cor.vec.cis) <- c('cor.CIS', 'pval.CIS')
cor.vec.cis$MetaProgram <- rownames(cor.vec.cis)

# Merge
cor.vec <- merge(cor.vec.dabtram, cor.vec.cocl2, by = 'MetaProgram')
cor.vec <- merge(cor.vec, cor.vec.cis, by = 'MetaProgram')
rownames(cor.vec) <- cor.vec$MetaProgram


mps <- c('Cell.Cycle...G1.S.rho', 'Cell.Cycle...G2.M.rho', 'Cell.Cycle.HMG.rich.rho', 'Chromatin.rho',
         'EMT.I.rho', 'EMT.II.rho', 'EMT.III.rho', 'EMT.IV.rho', 'Epithelial.Senescence.rho', 
         'Hypoxia.rho', 'Interferon.MHC.II..I..rho', 'Interferon.MHC.II..II..rho', 'MYC.rho', 
         'Proteasomal.degradation.rho', 'Protein.maturation.rho', 'Respiration.rho', 'Secreted.I.rho',
         'Secreted.II.rho', 'Skin.pigmentation.rho', 'Stress..in.vitro..rho', 'Stress.rho',
         'Translation.initiation.rho', 'Unfolded.protein.response.rho')

# =============================================================================
# Plot
# =============================================================================
# pheatmap::pheatmap(cor.vec[mps, c('cor.DABTRAM', 'cor.COCL2', 'cor.CIS')], 
#                    breaks = seq(-0.8, 0.8, length.out = 101),
#                    fontsize = 10, fontsize_row = 10, fontsize_col = 10,
#                    show_colnames = T, show_rownames = T,
#                    cluster_rows = F, cluster_cols = F, 
#                    border_color = 'black', 
#                    legend = F,
#                    cellwidth = 20, cellheight = 20,
#                    filename = paste0(figure_dir, 'Gavish_MP_correlation_FB_heatmap.pdf'))

pdf(paste0(figure_dir, 'Supp_Gavish_MP_correlation_FB_heatmap.pdf'), width = 4, height = 8)
pheatmap(cor.vec[mps, c('cor.DABTRAM', 'cor.COCL2', 'cor.CIS')], 
         breaks = seq(-0.8, 0.8, length.out = 101),
         fontsize = 10, fontsize_row = 10, fontsize_col = 10,
         show_colnames = T, show_rownames = T,
         cluster_rows = F, cluster_cols = F, 
         border_color = 'black', 
         cellwidth = 20, cellheight = 20)
dev.off()



df.bias.DABTRAM.1 <- read.csv(paste0(out_dir2, 'adapting_bias_thres_0_DABTRAM.csv'))
df.bias.COCL2.1 <- read.csv(paste0(out_dir2, 'adapting_bias_thres_0_COCL2.csv'))
df.bias.CIS.1 <- read.csv(paste0(out_dir2, 'adapting_bias_thres_0_CIS.csv'))

df.bias.DABTRAM.1 <- df.bias.DABTRAM.1[, c('cell_id', 'bias')]
df.bias.COCL2.1 <- df.bias.COCL2.1[, c('cell_id', 'bias')]
df.bias.CIS.1 <- df.bias.CIS.1[, c('cell_id', 'bias')]

colnames(df.bias.DABTRAM.1) <- c('cell_id', 'bias.DABTRAM')
colnames(df.bias.COCL2.1) <- c('cell_id', 'bias.COCL2')
colnames(df.bias.CIS.1) <- c('cell_id', 'bias.CIS')

df.bias.all <- merge(df.bias.DABTRAM.1, df.bias.COCL2.1, by = 'cell_id')
df.bias.all <- merge(df.bias.all, df.bias.CIS.1, by = 'cell_id')

pheatmap::pheatmap(df.bias.all[, c('bias.DABTRAM', 'bias.COCL2', 'bias.CIS')], 
                  fontsize = 10, 
                  fontsize_row = 10, 
                  fontsize_col = 10,
                  show_colnames = T,
                  show_rownames = F,
                  cluster_rows = T,
                  cluster_cols = F)

ggpairs(df.bias.all[, c('bias.DABTRAM', 'bias.COCL2', 'bias.CIS')], method = 'spearman')

        