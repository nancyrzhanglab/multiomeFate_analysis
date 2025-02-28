rm(list = ls())
library(tidyverse)
library(ggplot2)
library(GGally)
library(hdrcde)
library(ggdensity)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(gridExtra)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

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
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/saver_cor_vec.RData')

dabtram_d0_saver_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 
cocl2_d0_saver_cor_vec <- saver_cor_vec[['cocl2_d0_saver_cor_vec']] 
cis_d0_saver_cor_vec <- saver_cor_vec[['cis_d0_saver_cor_vec']] 

dabtram_d10_saver_cor_vec <- saver_cor_vec[['dabtram_d10_saver_cor_vec']] 
cocl2_d10_saver_cor_vec <- saver_cor_vec[['cocl2_d10_saver_cor_vec']] 
cis_d10_saver_cor_vec <- saver_cor_vec[['cis_d10_saver_cor_vec']] 

dabtram_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/dabtram_saver_cor_SMS.csv')
cocl2_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cocl2_saver_cor_SMS.csv')
cis_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cis_saver_cor_SMS.csv')

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_saver_cor_vec) <- paste0(colnames(dabtram_d0_saver_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_saver_cor_vec) <- paste0(colnames(cocl2_d0_saver_cor_vec), '.COCL2_d0')
colnames(cis_d0_saver_cor_vec) <- paste0(colnames(cis_d0_saver_cor_vec), '.CIS_d0')

colnames(dabtram_d10_saver_cor_vec) <- paste0(colnames(dabtram_d10_saver_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_saver_cor_vec) <- paste0(colnames(cocl2_d10_saver_cor_vec), '.COCL2_d10')
colnames(cis_d10_saver_cor_vec) <- paste0(colnames(cis_d10_saver_cor_vec), '.CIS_d10')

# Day0
d0_cor <- merge(dabtram_d0_saver_cor_vec, cocl2_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

d0_cor <- merge(d0_cor, cis_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

# Day10
d10_cor <- merge(dabtram_d10_saver_cor_vec, cocl2_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

d10_cor <- merge(d10_cor, cis_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

# Filter for genes of interest
dabtram_annotations <- dabtram_annotations[dabtram_annotations$Gene.of.interest == 'Yes', ]
cocl2_annotations <- cocl2_annotations[cocl2_annotations$Gene.of.interest == 'Yes', ]
cis_annotations <- cis_annotations[cis_annotations$Gene.of.interest == 'Yes', ]

d0_cor_anno <- d0_cor[rownames(d0_cor) %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]

# ==============================================================================
# Plotting
# ==============================================================================

ggpairs(d0_cor[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0')], mapping = aes(alpha = 0.7, color = 'pink')) +
  theme_Publication()

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_RNA_saver_FP_day0.png', 
       width = 8, height = 8, dpi = 300)

ggpairs(d10_cor[, c('correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')], mapping = aes(alpha = 0.7, color = 'pink')) +
  theme_Publication()

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_RNA_saver_FP_day10.png', 
       width = 8, height = 8, dpi = 300)

# ==============================================================================
# Plotting individual ones
# ==============================================================================

ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.COCL2_d0)) +
  geom_point(color = 'pink') +
  theme_Publication()

d0_cor <- d0_cor %>% drop_na()

# scatter plot in density
ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.COCL2_d0)) +
  geom_point(color = 'pink') +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) +
  # geom_density_2d(color = 'pink') +
  expand_limits(x = range(d0_cor$correlation.DABTRAM_d0), y = range(d0_cor$correlation.COCL2_d0)) +
  theme_Publication()

p1 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p2 <- ggplot(d0_cor, aes(x = correlation.CIS_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()


p3 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.CIS_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p4 <- grid.arrange(p1, p2, p3, ncol = 3)
ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_RNA_saver_FP_day0_panel.pdf', 
       p4, 
       width = 8, height = 3, dpi = 300)


d10_cor <- d10_cor %>% drop_na()

p4 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p5 <- ggplot(d10_cor, aes(x = correlation.CIS_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()


p6 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10, y = correlation.CIS_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p7 <- grid.arrange(p4, p5, p6, ncol = 3)

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/correlation_RNA_saver_FP_day10_panel.pdf', 
       p7, 
       width = 8, height = 3, dpi = 300)

# =============================================================================
# Plot
# =============================================================================
d0_cor$gene <- rownames(d0_cor)
d10_cor$gene <- rownames(d10_cor)

dabtram_genes <- 
  union(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3"),
            c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44"))

dabtram_genes <- c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")


cocl2_genes <- sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
               "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
               "CAV1"))
cis_genes <- sort(c("YY1AP1", "LGALS3", "MCF2L", "TIMM50", "AC207130.1",
             "SLC25A6", "EIF3L", "CTSD", "NQO1", "HNMT", "ZFYVE16",
             "PHACTR1", "TNFRSF14", "RAI14", "TRPM1", "HIST1H1C",
             "HIST2H2AC", "SPARC", "TRIM63", "TUBA1B", "HIST1H1A",
             "HIST1H1D", "PYCARD", "FSTL1", "DCT", "CTSK", "HIST1H4C",
             "GDF15", "HIST1H1B"))

keygenes.dabtram.d0 <- d0_cor[d0_cor$gene %in% dabtram_genes, 'correlation.DABTRAM_d0']
keygenes.cocl2.d0 <- d0_cor[d0_cor$gene %in% cocl2_genes, 'correlation.COCL2_d0']
keygenes.cis.d0 <- d0_cor[d0_cor$gene %in% cis_genes, 'correlation.CIS_d0']

keygenes.dabtram.d10 <- d10_cor[d10_cor$gene %in% dabtram_genes, 'correlation.DABTRAM_d10']
keygenes.cocl2.d10 <- d10_cor[d10_cor$gene %in% cocl2_genes, 'correlation.COCL2_d10']
keygenes.cis.d10 <- d10_cor[d10_cor$gene %in% cis_genes, 'correlation.CIS_d10']

x_min <- min(c(d0_cor$correlation.DABTRAM_d0, d0_cor$correlation.COCL2_d0, 
               d0_cor$correlation.CIS_d0, d10_cor$correlation.DABTRAM_d10, 
               d10_cor$correlation.COCL2_d10, d10_cor$correlation.CIS_d10), rm.na = T)
x_min <- -1
x_max <- max(c(d0_cor$correlation.DABTRAM_d0, d0_cor$correlation.COCL2_d0, 
               d0_cor$correlation.CIS_d0, d10_cor$correlation.DABTRAM_d10, 
               d10_cor$correlation.COCL2_d10, d10_cor$correlation.CIS_d10), rm.na = T)
x_max <- 1

p1 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.dabtram.d0, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, gene %in% dabtram_genes), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(x_min, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (DABTRAM)')
p1

p2 <- ggplot(d0_cor, aes(x = correlation.COCL2_d0)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.cocl2.d0, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, gene %in% cocl2_genes), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (COCL2)')

p3 <- ggplot(d0_cor, aes(x = correlation.CIS_d0)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.cis.d0, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, gene %in% cis_genes), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (CIS)')

p4 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.dabtram.d10, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, gene %in% dabtram_genes), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (DABTRAM)')

p5 <- ggplot(d10_cor, aes(x = correlation.COCL2_d10)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.cocl2.d10, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, gene %in% cocl2_genes), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (COCL2)')

p6 <- ggplot(d10_cor, aes(x = correlation.CIS_d10)) +
  geom_density(fill = 'grey') +
  geom_vline(xintercept = keygenes.cis.d10, color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, gene %in% cis_genes), 
                           aes(label = gene, y = 3)) +
  theme_Publication() +
  xlim(-1, x_max) +
  labs(ylab = 'Density', xlab = 'Correlation (CIS)')

p7 <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 1)

ggsave('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/geneSaver_correlation_density.pdf', 
       plot = p7, width = 6, height = 8)

top <- d10_cor %>% 
  mutate(padjust = p.adjust(p.value.DABTRAM_d10, method = 'BH')) %>%
  filter(correlation.DABTRAM_d10 > 0) %>% 
  filter(padjust < 0.05) %>%
  arrange(correlation.DABTRAM_d10) %>%
  head(1) %>% 
  pull(correlation.DABTRAM_d10)

bottom <- d10_cor %>% 
  mutate(padjust = p.adjust(p.value.DABTRAM_d10, method = 'BH')) %>%
  filter(correlation.DABTRAM_d10 < 0) %>% 
  filter(padjust < 0.05) %>%
  arrange(correlation.DABTRAM_d10) %>%
  tail(1) %>% 
  pull(correlation.DABTRAM_d10)

d10_cor <- d10_cor[order(d10_cor$correlation.DABTRAM_d10, decreasing = T),]
d10_cor$order <- 1:nrow(d10_cor)
ggplot(d10_cor, aes(x = order, y = correlation.DABTRAM_d10)) +
  geom_point() +
  geom_point(data = subset(d10_cor, gene %in% dabtram_genes), color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, gene %in% dabtram_genes), color = 'red',
                           aes(label = gene), box.padding = unit(0.3, 'lines'),
                           point.padding = unit(0.5, 'lines'),
                           max.overlaps = 50) +
  geom_hline(yintercept = top, color = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = bottom, color = 'gray', linetype = 'dashed') +
  ggtitle('GEX') +
  theme_Publication()



