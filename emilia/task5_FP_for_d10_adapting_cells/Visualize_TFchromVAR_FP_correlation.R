rm(list = ls())

library(ggplot2)
library(GGally)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(gridExtra)

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
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
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

# =============================================================================
# Read data
# =============================================================================
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

load(paste0(out_dir, 'TFchromVAR_on_day0_cor_vec_DABTRAM.RData'))
cor_vec.DABTRAM <- chromVAR_cor_vec

load(paste0(out_dir, 'TFchromVAR_on_day0_cor_vec_COCL2.RData'))
cor_vec.COCL2 <- chromVAR_cor_vec

load(paste0(out_dir, 'TFchromVAR_on_day0_cor_vec_CIS.RData'))
cor_vec.CIS <- chromVAR_cor_vec

# =============================================================================
# Wrangle
# =============================================================================
cor_vec.DABTRAM$TF <- rownames(cor_vec.DABTRAM)
cor_vec.COCL2$TF <- rownames(cor_vec.COCL2)
cor_vec.CIS$TF <- rownames(cor_vec.CIS)

colnames(cor_vec.DABTRAM) <- c('cor.DATBRAM', 'p_val.DABTRAM', 'TF')
colnames(cor_vec.COCL2) <- c('cor.COCL2', 'p_val.COCL2', 'TF')
colnames(cor_vec.CIS) <- c('cor.CIS', 'p_val.CIS', 'TF')

comp_df <- merge(cor_vec.DABTRAM, cor_vec.COCL2, by = 'TF')
comp_df <- merge(comp_df, cor_vec.CIS, by = 'TF')

ggpairs(comp_df[, c('cor.DATBRAM', 'cor.COCL2', 'cor.CIS')], 
        lower = list(continuous = 'points'))

tfs <- comp_df$TF
tfs_toplot <- tfs[grepl('JUN', tfs) | grepl('FOS', tfs) | grepl('SOX10', tfs) | grepl('MITF', tfs) | grepl('TEAD', tfs) ]
tfs_toplot <- tfs_toplot[!grepl('(var.2)', tfs_toplot)]

# =============================================================================
# Plot pair-wise comparisons
# =============================================================================
x_min <- min(comp_df[, c('cor.DATBRAM', 'cor.COCL2', 'cor.CIS')]) + 0.08
x_max <- max(comp_df[, c('cor.DATBRAM', 'cor.COCL2', 'cor.CIS')]) + 0.08

p1 <- ggplot(comp_df, aes(x = cor.DATBRAM, y = cor.COCL2)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.25, label.y = 0.45) +
  xlim(c(x_min, x_max)) +
  ylim(c(x_min, x_max)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p2 <- ggplot(comp_df, aes(x = cor.CIS, y = cor.COCL2)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.25, label.y = 0.45) +
  xlim(c(x_min, x_max)) +
  ylim(c(x_min, x_max)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p3 <- ggplot(comp_df, aes(x = cor.DATBRAM, y = cor.CIS)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.25, label.y = 0.45) +
  xlim(c(x_min, x_max)) +
  ylim(c(x_min, x_max)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p4 <- grid.arrange(p1, p2, p3, ncol = 3)
ggsave(paste0(out_dir, 'TFchromVAR_correlation_day0_panel.pdf'), 
       p4, 
       width = 8, height = 3, dpi = 300)

# =============================================================================
# Plot
# =============================================================================
keyTFs.dabtram <- comp_df[comp_df$TF %in% tfs_toplot, 'cor.DATBRAM']
keyTFs.cocl2 <- comp_df[comp_df$TF %in% tfs_toplot, 'cor.COCL2']
keyTFs.cis <- comp_df[comp_df$TF %in% tfs_toplot, 'cor.CIS']

x_min <- min(c(comp_df$cor.DATBRAM, comp_df$cor.COCL2, comp_df$cor.CIS))
x_max <- max(c(comp_df$cor.DATBRAM, comp_df$cor.COCL2, comp_df$cor.CIS))

# DABTRAM
p1 <- ggplot(comp_df, aes(x = cor.DATBRAM)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.dabtram, color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(x_min, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (DABTRAM)')

# COCL2
p2 <- ggplot(comp_df, aes(x = cor.COCL2)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.cocl2, color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(x_min, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (COCL2)')

# CIS
p3 <- ggplot(comp_df, aes(x = cor.CIS)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keyTFs.cis, color = 'red') +
  ggrepel::geom_text_repel(data = subset(comp_df, TF %in% tfs_toplot), 
                           aes(label = TF, y = 3)) +
  theme_Publication() +
  xlim(x_min, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (CIS)')

p4 <- grid.arrange(p1, p2, p3, ncol = 1)

ggsave(paste0(out_dir, 'TFchromVAR_correlation_density.pdf'), plot = p4, width = 6, height = 4)
