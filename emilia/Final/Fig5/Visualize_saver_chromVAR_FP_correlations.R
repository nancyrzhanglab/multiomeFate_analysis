rm(list = ls())
library(tidyverse)
library(data.table)
library(ggplot2)
library(GGally)
library(hdrcde)
library(ggdensity)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggpubr)
library(ggrepel)
library(gridExtra)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

theme_Publication<- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = base_size),
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

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig5/'


##############################################################################################
## Fig 5B: Density plot of saver FP correlations
##############################################################################################

# ==============================================================================
# Read data
# ==============================================================================

load(paste0(result_dir, 'saver_cor_vec.RData'))

dabtram_d0_saver_cor_vec <- saver_cor_vec[['dabtram_d0_saver_cor_vec']] 
cocl2_d0_saver_cor_vec <- saver_cor_vec[['cocl2_d0_saver_cor_vec']] 
cis_d0_saver_cor_vec <- saver_cor_vec[['cis_d0_saver_cor_vec']] 

dabtram_d10_saver_cor_vec <- saver_cor_vec[['dabtram_d10_saver_cor_vec']] 
cocl2_d10_saver_cor_vec <- saver_cor_vec[['cocl2_d10_saver_cor_vec']] 
cis_d10_saver_cor_vec <- saver_cor_vec[['cis_d10_saver_cor_vec']] 

# ==============================================================================
# Format data
# ==============================================================================
colnames(dabtram_d0_saver_cor_vec) <- paste0(colnames(dabtram_d0_saver_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_saver_cor_vec) <- paste0(colnames(cocl2_d0_saver_cor_vec), '.COCL2_d0')
colnames(cis_d0_saver_cor_vec) <- paste0(colnames(cis_d0_saver_cor_vec), '.CIS_d0')

colnames(dabtram_d10_saver_cor_vec) <- paste0(colnames(dabtram_d10_saver_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_saver_cor_vec) <- paste0(colnames(cocl2_d10_saver_cor_vec), '.COCL2_d10')
colnames(cis_d10_saver_cor_vec) <- paste0(colnames(cis_d10_saver_cor_vec), '.CIS_d10')

# make a master dataframe for Day0
d0_cor.RNA <- merge(dabtram_d0_saver_cor_vec, cocl2_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor.RNA) <- d0_cor.RNA$Row.names
d0_cor.RNA <- d0_cor.RNA |> select(-Row.names)

d0_cor.RNA <- merge(d0_cor.RNA, cis_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor.RNA) <- d0_cor.RNA$Row.names
d0_cor.RNA <- d0_cor.RNA |> select(-Row.names)

d0_cor.RNA$gene <- rownames(d0_cor.RNA)

# make a master dataframe for Day10
d10_cor.RNA <- merge(dabtram_d10_saver_cor_vec, cocl2_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor.RNA) <- d10_cor.RNA$Row.names
d10_cor.RNA <- d10_cor.RNA |> select(-Row.names)

d10_cor.RNA <- merge(d10_cor.RNA, cis_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor.RNA) <- d10_cor.RNA$Row.names
d10_cor.RNA <- d10_cor.RNA |> select(-Row.names)

d10_cor.RNA$gene <- rownames(d10_cor.RNA)

# =============================================================================
# Plot density
# =============================================================================

# Plot params
x_min <- -1
x_max <- 1


# Define key genes
res_genes <- c("WNT5A", "AXL", "EGFR", "JUN", "NGFR", "PCNA")

# Pick out res genes
keygenes.dabtram.d0 <- d0_cor.RNA[d0_cor.RNA$gene %in% res_genes, 'correlation.DABTRAM_d0']
keygenes.cocl2.d0 <- d0_cor.RNA[d0_cor.RNA$gene %in% res_genes, 'correlation.COCL2_d0']
keygenes.cis.d0 <- d0_cor.RNA[d0_cor.RNA$gene %in% res_genes, 'correlation.CIS_d0']

keygenes.dabtram.d10 <- d10_cor.RNA[d10_cor.RNA$gene %in% res_genes, 'correlation.DABTRAM_d10']
keygenes.cocl2.d10 <- d10_cor.RNA[d10_cor.RNA$gene %in% res_genes, 'correlation.COCL2_d10']
keygenes.cis.d10 <- d10_cor.RNA[d10_cor.RNA$gene %in% res_genes, 'correlation.CIS_d10']

# DABTRAM day0
d0_cor.RNA <- d0_cor.RNA %>% arrange(correlation.DABTRAM_d0)
d0_cor.RNA$order.DABTRAM <- 1:nrow(d0_cor.RNA)
p1 <- ggplot(d0_cor.RNA, aes(x = order.DABTRAM, y = correlation.DABTRAM_d0)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = d0_cor.RNA[d0_cor.RNA$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor.RNA, gene %in% res_genes), 
                            box.padding = 0.1, 
                            max.overlaps = Inf,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

# COCL2 day0
d0_cor.RNA <- d0_cor.RNA %>% arrange(correlation.COCL2_d0)
d0_cor.RNA$order.COCL2 <- 1:nrow(d0_cor.RNA)
p2 <- ggplot(d0_cor.RNA, aes(x = order.COCL2, y = correlation.COCL2_d0)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = d0_cor.RNA[d0_cor.RNA$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor.RNA, gene %in% res_genes), 
                           box.padding = 0.1, 
                           max.overlaps = Inf,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

# CIS day0
d0_cor.RNA <- d0_cor.RNA %>% arrange(correlation.CIS_d0)
d0_cor.RNA$order.CIS <- 1:nrow(d0_cor.RNA)
p3 <- ggplot(d0_cor.RNA, aes(x = order.CIS, y = correlation.CIS_d0)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = d0_cor.RNA[d0_cor.RNA$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor.RNA, gene %in% res_genes), 
                            box.padding = 0.1, 
                            max.overlaps = Inf,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

p4 <- grid.arrange(p1, p2, p3, ncol = 3)
ggsave(filename = paste0(figure_dir, 'SuppFig5A.D0.RNA.pdf'), p4, width = 8.5, height = 2.5)




# DABTRAM day10
d10_cor.RNA <- d10_cor.RNA %>% arrange(correlation.DABTRAM_d10)
d10_cor.RNA$order.DABTRAM <- 1:nrow(d10_cor.RNA)
p5 <- ggplot(d10_cor.RNA, aes(x = order.DABTRAM, y = correlation.DABTRAM_d10)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = d10_cor.RNA[d10_cor.RNA$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_label_repel(data = subset(d10_cor.RNA, gene %in% res_genes), 
                           box.padding = 0.5, 
                           max.overlaps = Inf,
                           nudge_x = -2,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

# COCL2 day10
d10_cor.RNA <- d10_cor.RNA %>% arrange(correlation.COCL2_d10)
d10_cor.RNA$order.COCL2 <- 1:nrow(d10_cor.RNA)
p6 <- ggplot(d10_cor.RNA, aes(x = order.COCL2, y = correlation.COCL2_d10)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = d10_cor.RNA[d10_cor.RNA$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_label_repel(data = subset(d10_cor.RNA, gene %in% res_genes), 
                           box.padding = 0.5, 
                           max.overlaps = Inf,
                           nudge_x = -10,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')


# CIS day10
d10_cor.RNA <- d10_cor.RNA %>% arrange(correlation.CIS_d10)
d10_cor.RNA$order.CIS <- 1:nrow(d10_cor.RNA)
p7 <- ggplot(d10_cor.RNA, aes(x = order.CIS, y = correlation.CIS_d10)) +
  geom_point(size = 1, color = '#FFDBDB') +
  geom_point(data = d10_cor.RNA[d10_cor.RNA$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_label_repel(data = subset(d10_cor.RNA, gene %in% res_genes), 
                           box.padding = 0.5,
                           max.overlaps = Inf,
                           nudge_x = -10,
                           aes(label = gene)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

p8 <- grid.arrange(p5, p6, p7, ncol = 3)

ggsave(filename = paste0(figure_dir, 'Fig5B.RNA.pdf'), p8, width = 8.5, height = 2.5)

##############################################################################################
## Fig 5D: Comparison of correlations
##############################################################################################

# Day0 correlations
p9 <- ggplot(d0_cor.RNA, aes(x = correlation.DABTRAM_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication() +
  theme(legend.position = "none")

p10 <- ggplot(d0_cor.RNA, aes(x = correlation.CIS_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication() +
  theme(legend.position = "none")


p11 <- ggplot(d0_cor.RNA, aes(x = correlation.DABTRAM_d0, y = correlation.CIS_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication() +
  theme(legend.position = "none")

p12 <- grid.arrange(p9, p10, p11, ncol = 3)

ggsave(filename = paste0(figure_dir, 'Fig5D.D0.RNA.pdf'), p12, width = 9, height = 2.8)

# D10 correlations
p13 <- ggplot(d10_cor.RNA, aes(x = correlation.DABTRAM_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p14 <- ggplot(d10_cor.RNA, aes(x = correlation.CIS_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()


p15 <- ggplot(d10_cor.RNA, aes(x = correlation.DABTRAM_d10, y = correlation.CIS_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication() +

p16 <- grid.arrange(p13, p14, p15, ncol = 1)


ggsave(filename = paste0(figure_dir, 'Fig5D.D10.RNA.pdf'), p16, width = 3, height = 9)

p13.cl <- p13 + 
  xlab('') +
  ylab('') +
  theme(legend.position = "none")

p14.cl <- p14 + 
  xlab('') +
  ylab('') +
  theme(legend.position = "none")

p15.cl <- p15 +
  xlab('') +
  ylab('') +
  theme(legend.position = "none")

p16.cl <- grid.arrange(p13.cl, p14.cl, p15.cl, ncol = 1)
ggsave(filename = paste0(figure_dir, 'Fig5D.D10.RNA.clean.pdf'), p16.cl, width = 2.8, height = 7.5)

##############################################################################################
## Supp: Sample gene vs FP scatter plots
##############################################################################################
# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

remove_unassigned_cells <- TRUE

all_data[['saver']] <- all_data_saver
all_data@misc <- all_data_fatepotential

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}
# Define key genes
res_genes <- c("WNT5A", "AXL", "EGFR", "JUN", "NGFR", "PCNA", "RPL8")
# =============================================================================
# Wrangle data
# =============================================================================

metadat.all_data <- all_data@meta.data
metadat.all_data$cell_id <- rownames(metadat.all_data)

saver.all_data <- all_data[['saver']]@scale.data

saver.all_data.toplot <- saver.all_data[res_genes, ]
saver.all_data.toplot <- t(saver.all_data.toplot)
saver.all_data.toplot <- as.data.frame(saver.all_data.toplot)
saver.all_data.toplot$cell_id <- rownames(saver.all_data.toplot)

# fatepotential.all_data <- metadat.all_data[, c('cell_id', 'dataset',
#                                                'fatepotential_CIS_d0_d10', 'fatepotential_CIS_d10_w5', 
#                                                'fatepotential_COCL2_d0_d10', 'fatepotential_COCL2_d10_w5',
#                                                'fatepotential_DABTRAM_d0_d10', 'fatepotential_DABTRAM_d10_w5')]
fatepotential.all_data <- metadat.all_data[, c('cell_id', 'dataset', 'fatepotential_CIS_d10_w5', 
                                               'fatepotential_COCL2_d10_w5', 'fatepotential_DABTRAM_d10_w5')]

cor.to.plot <- merge(fatepotential.all_data, saver.all_data.toplot, by = 'cell_id')

# =============================================================================
# Plot data
# =============================================================================
# cor.to.plot.melt <- melt(cor.to.plot, id.vars = c('cell_id', 'dataset', 'fatepotential_CIS_d0_d10', 'fatepotential_CIS_d10_w5', 
#                                                   'fatepotential_COCL2_d0_d10', 'fatepotential_COCL2_d10_w5',
#                                                   'fatepotential_DABTRAM_d0_d10', 'fatepotential_DABTRAM_d10_w5'))
# colnames(cor.to.plot.melt)[9] <- 'Gene'
# colnames(cor.to.plot.melt)[10] <- 'Saver'
# cor.to.plot.melt <- melt(cor.to.plot.melt, id.vars = c('cell_id', 'dataset', 'Gene', 'Saver'))
# colnames(cor.to.plot.melt)[5] <- 'FatePotential.index'
# colnames(cor.to.plot.melt)[6] <- 'FatePotential.estimate'

cor.to.plot.melt <- melt(cor.to.plot, id.vars = c('cell_id', 'dataset', 'fatepotential_CIS_d10_w5',
                                                  'fatepotential_COCL2_d10_w5','fatepotential_DABTRAM_d10_w5'))
colnames(cor.to.plot.melt)[6] <- 'Gene'
colnames(cor.to.plot.melt)[7] <- 'Saver'
cor.to.plot.melt <- melt(cor.to.plot.melt, id.vars = c('cell_id', 'dataset', 'Gene', 'Saver'))
colnames(cor.to.plot.melt)[5] <- 'FatePotential.index'
colnames(cor.to.plot.melt)[6] <- 'FatePotential.estimate'

cor.to.plot.melt <- cor.to.plot.melt %>% drop_na()
cor.to.plot.melt$FatePotential.index <- factor(cor.to.plot.melt$FatePotential.index, levels = c('fatepotential_DABTRAM_d10_w5',
                                                                                                'fatepotential_COCL2_d10_w5', 
                                                                                                'fatepotential_CIS_d10_w5'))
cor.to.plot.melt.ctrl <- cor.to.plot.melt[cor.to.plot.melt$Gene == 'RPL8', ]
cor.to.plot.melt <- cor.to.plot.melt[cor.to.plot.melt$Gene != 'RPL8', ]

cor.to.plot.melt$FatePotential.index <- gsub('fatepotential_', '', cor.to.plot.melt$FatePotential.index)
cor.to.plot.melt$FatePotential.index <- gsub('_d10_w5', '', cor.to.plot.melt$FatePotential.index)
cor.to.plot.melt$FatePotential.index <- factor(cor.to.plot.melt$FatePotential.index, levels = c('DABTRAM',
                                                                                                'COCL2', 
                                                                                                'CIS'))

p1 <- ggplot(cor.to.plot.melt, aes(x = Saver, y = FatePotential.estimate)) +
  geom_point(color = 'gray', size = 1, alpha = 0.5) +
  # geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.sep = "\n", na.rm = T, label.x = 0, label.y = 0) +
  xlab('GEX (Saver)') +
  ylab('D10-to-W5 Growth Potential') +
  facet_grid(FatePotential.index ~ Gene, scales = 'free') +
  coord_cartesian(xlim = c(-1, 5)) +
  theme_bw()

cor.to.plot.melt.ctrl$FatePotential.index <- gsub('fatepotential_', '', cor.to.plot.melt.ctrl$FatePotential.index)
cor.to.plot.melt.ctrl$FatePotential.index <- gsub('_d10_w5', '', cor.to.plot.melt.ctrl$FatePotential.index)
cor.to.plot.melt.ctrl$FatePotential.index <- factor(cor.to.plot.melt.ctrl$FatePotential.index, levels = c('DABTRAM',
                                                                                                'COCL2', 
                                                                                                'CIS'))

p2 <- ggplot(cor.to.plot.melt.ctrl, aes(x = Saver, y = FatePotential.estimate)) +
  geom_point(color = 'gray', size = 1, alpha = 0.5) +
  # geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.sep = "\n", na.rm = T, label.x = -7.5, label.y = 0) +
  xlab('GEX (Saver)') +
  ylab('D10-to-W5 Growth Potential') +
  facet_grid(FatePotential.index ~ Gene, scales = 'free') +
  # coord_cartesian(xlim = c(-1, 5)) +
  theme_bw()

p3 <- ggarrange(p1, p2, widths = c(6, 1.5))

ggsave(p3, filename = paste0(figure_dir, 'Supp_Geve_vs_FatePotential.png'), 
       width = 9.5, height = 4, dpi = 600)

##############################################################################################
## Fig 5B: Density plot of chromVAR FP correlations
##############################################################################################
# ==============================================================================
# Read data
# ==============================================================================
load('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/chromVAR_cor_vec.RData')

dabtram_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d0_chromVAR_cor_vec']] 
cocl2_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d0_chromVAR_cor_vec']] 
cis_d0_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d0_chromVAR_cor_vec']] 

dabtram_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['dabtram_d10_chromVAR_cor_vec']] 
cocl2_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cocl2_d10_chromVAR_cor_vec']] 
cis_d10_chromVAR_cor_vec <- chromVAR_cor_vec[['cis_d10_chromVAR_cor_vec']] 

# ==============================================================================
# Wrangle data
# ==============================================================================
colnames(dabtram_d0_chromVAR_cor_vec) <- paste0(colnames(dabtram_d0_chromVAR_cor_vec), '.DABTRAM_d0')
colnames(cocl2_d0_chromVAR_cor_vec) <- paste0(colnames(cocl2_d0_chromVAR_cor_vec), '.COCL2_d0')
colnames(cis_d0_chromVAR_cor_vec) <- paste0(colnames(cis_d0_chromVAR_cor_vec), '.CIS_d0')

colnames(dabtram_d10_chromVAR_cor_vec) <- paste0(colnames(dabtram_d10_chromVAR_cor_vec), '.DABTRAM_d10')
colnames(cocl2_d10_chromVAR_cor_vec) <- paste0(colnames(cocl2_d10_chromVAR_cor_vec), '.COCL2_d10')
colnames(cis_d10_chromVAR_cor_vec) <- paste0(colnames(cis_d10_chromVAR_cor_vec), '.CIS_d10')

# Day0
d0_cor.chromVAR <- merge(dabtram_d0_chromVAR_cor_vec, cocl2_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor.chromVAR) <- d0_cor.chromVAR$Row.names
d0_cor.chromVAR <- d0_cor.chromVAR |> dplyr::select(-`Row.names`)

d0_cor.chromVAR <- merge(d0_cor.chromVAR, cis_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor.chromVAR) <- d0_cor.chromVAR$Row.names
d0_cor.chromVAR <- d0_cor.chromVAR |> dplyr::select(-`Row.names`)


# Day10
d10_cor.chromVAR <- merge(dabtram_d10_chromVAR_cor_vec, cocl2_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor.chromVAR) <- d10_cor.chromVAR$Row.names
d10_cor.chromVAR <- d10_cor.chromVAR |> dplyr::select(-Row.names)

d10_cor.chromVAR <- merge(d10_cor.chromVAR, cis_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor.chromVAR) <- d10_cor.chromVAR$Row.names
d10_cor.chromVAR <- d10_cor.chromVAR |> dplyr::select(-Row.names)


d0_cor.chromVAR$TF <- rownames(d0_cor.chromVAR)
d10_cor.chromVAR$TF <- rownames(d10_cor.chromVAR)
tfs <- d0_cor.chromVAR$TF
tfs_toplot <- tfs[grepl('JUN', tfs) | grepl('FOS', tfs) | grepl('SOX10', tfs) | grepl('MITF', tfs) | grepl('TEAD', tfs) | grepl('STAT', tfs) | grepl('IRF3', tfs) ]
tfs_toplot <- tfs_toplot[!grepl('(var.2)', tfs_toplot)]
tfs_toplot <- c('JUNB', 'FOS', 'FOSL2', 'TEAD2', 'TEAD1', 'JUN', 'STAT1', 'MITF', 'SOX10')

keyTFs.dabtram.d0 <- d0_cor.chromVAR[d0_cor.chromVAR$TF %in% tfs_toplot, 'correlation.DABTRAM_d0']
keyTFs.cocl2.d0 <- d0_cor.chromVAR[d0_cor.chromVAR$TF %in% tfs_toplot, 'correlation.COCL2_d0']
keyTFs.cis.d0 <- d0_cor.chromVAR[d0_cor.chromVAR$TF %in% tfs_toplot, 'correlation.CIS_d0']

keyTFs.dabtram.d10 <- d10_cor.chromVAR[d10_cor.chromVAR$TF %in% tfs_toplot, 'correlation.DABTRAM_d10']
keyTFs.cocl2.d10 <- d10_cor.chromVAR[d10_cor.chromVAR$TF %in% tfs_toplot, 'correlation.COCL2_d10']
keyTFs.cis.d10 <- d10_cor.chromVAR[d10_cor.chromVAR$TF %in% tfs_toplot, 'correlation.CIS_d10']

x_min <- -1
x_max <- 1

# DABTRAM day0
d0_cor.chromVAR <- d0_cor.chromVAR %>% arrange(correlation.DABTRAM_d0)
d0_cor.chromVAR$order.DABTRAM <- 1:nrow(d0_cor.chromVAR)
p1.TF <- ggplot(d0_cor.chromVAR, aes(x = order.DABTRAM, y = correlation.DABTRAM_d0)) +
  geom_point(size = 1, color = '#B5DFB7') +
  geom_point(data = d0_cor.chromVAR[d0_cor.chromVAR$TF %in% tfs_toplot, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor.chromVAR, TF %in% tfs_toplot), 
                           aes(label = TF),
                           max.overlaps = Inf) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'TF rank')


# COCL2 day0
d0_cor.chromVAR <- d0_cor.chromVAR %>% arrange(correlation.COCL2_d0)
d0_cor.chromVAR$order.COCL2 <- 1:nrow(d0_cor.chromVAR)
p2.TF <- ggplot(d0_cor.chromVAR, aes(x = order.COCL2, y = correlation.COCL2_d0)) +
  geom_point(size = 1, color = '#B5DFB7') +
  geom_point(data = d0_cor.chromVAR[d0_cor.chromVAR$TF %in% tfs_toplot, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor.chromVAR, TF %in% tfs_toplot), 
                           aes(label = TF),
                           max.overlaps = Inf) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'TF rank')

# CIS day0
d0_cor.chromVAR <- d0_cor.chromVAR %>% arrange(correlation.CIS_d0)
d0_cor.chromVAR$order.CIS <- 1:nrow(d0_cor.chromVAR)
p3.TF <- ggplot(d0_cor.chromVAR, aes(x = order.CIS, y = correlation.CIS_d0)) +
  geom_point(size = 1, color = '#B5DFB7') +
  geom_point(data = d0_cor.chromVAR[d0_cor.chromVAR$TF %in% tfs_toplot, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor.chromVAR, TF %in% tfs_toplot), 
                           aes(label = TF),
                           max.overlaps = Inf) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'TF rank')

p4.TF <- grid.arrange(p1.TF, p2.TF, p3.TF, ncol = 3)
ggsave(filename = paste0(figure_dir, 'SuppFig5B.D0.TF.pdf'), p4.TF, width = 8.5, height = 2.5)


# DABTRAM day10
d10_cor.chromVAR <- d10_cor.chromVAR %>% arrange(correlation.DABTRAM_d10)
d10_cor.chromVAR$order.DABTRAM <- 1:nrow(d10_cor.chromVAR)
p5.TF <- ggplot(d10_cor.chromVAR, aes(x = order.DABTRAM, y = correlation.DABTRAM_d10)) +
  geom_point(size = 1, color = '#B5DFB7') +
  geom_point(data = d10_cor.chromVAR[d10_cor.chromVAR$TF %in% tfs_toplot, ], color = 'red') +
  ggrepel::geom_label_repel(data = subset(d10_cor.chromVAR, TF %in% tfs_toplot), 
                            box.padding = 0.5,
                            max.overlaps = Inf,
                           nudge_x = -2,
                           aes(label = TF)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'TF rank')

# COCL2 day10
d10_cor.chromVAR <- d10_cor.chromVAR %>% arrange(correlation.COCL2_d10)
d10_cor.chromVAR$order.COCL2 <- 1:nrow(d10_cor.chromVAR)
p6.TF <- ggplot(d10_cor.chromVAR, aes(x = order.COCL2, y = correlation.COCL2_d10)) +
  geom_point(size = 1, color = '#B5DFB7') +
  geom_point(data = d10_cor.chromVAR[d10_cor.chromVAR$TF %in% tfs_toplot, ], color = 'red') +
  ggrepel::geom_label_repel(data = subset(d10_cor.chromVAR, TF %in% tfs_toplot), 
                           nudge_x = -10,
                           box.padding = 0.5,
                           max.overlaps = Inf,
                           aes(label = TF)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'TF rank')


# CIS day10
d10_cor.chromVAR <- d10_cor.chromVAR %>% arrange(correlation.CIS_d10)
d10_cor.chromVAR$order.CIS <- 1:nrow(d10_cor.chromVAR)
p7.TF <- ggplot(d10_cor.chromVAR, aes(x = order.CIS, y = correlation.CIS_d10)) +
  geom_point(size = 1, color = '#B5DFB7') +
  geom_point(data = d10_cor.chromVAR[d10_cor.chromVAR$TF %in% tfs_toplot, ], color = 'red') +
  ggrepel::geom_label_repel(data = subset(d10_cor.chromVAR, TF %in% tfs_toplot), 
                           nudge_x = -10,
                           box.padding = 0.5,
                           max.overlaps = Inf,
                           aes(label = TF)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'TF rank')

p8.TF <- grid.arrange(p5.TF, p6.TF, p7.TF, ncol = 3)

ggsave(filename = paste0(figure_dir, 'Fig5B.TF.pdf'), p8.TF, width = 8.5, height = 2.5)

##############################################################################################
## Fig 5D: Comparison of correlations (TF)
##############################################################################################

# Day0 correlations
p9.TF <- ggplot(d0_cor.chromVAR, aes(x = correlation.DABTRAM_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication() +
  theme(legend.position = "none")

p10.TF <- ggplot(d0_cor.chromVAR, aes(x = correlation.CIS_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication() +
  theme(legend.position = "none")


p11.TF <- ggplot(d0_cor.chromVAR, aes(x = correlation.DABTRAM_d0, y = correlation.CIS_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "pearson", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication() +
  theme(legend.position = "none")

p12.TF <- grid.arrange(p9.TF, p10.TF, p11.TF, ncol = 3)

ggsave(filename = paste0(figure_dir, 'Fig5D.D0.TF.pdf'), p12.TF, width = 9, height = 2.8)


# D10 correlations
p13.TF <- ggplot(d10_cor.chromVAR, aes(x = correlation.DABTRAM_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p14.TF <- ggplot(d10_cor.chromVAR, aes(x = correlation.CIS_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.15, label.y = 0.8) +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()


p15.TF <- ggplot(d10_cor.chromVAR, aes(x = correlation.DABTRAM_d10, y = correlation.CIS_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.15) +
  xlim(c(-1, 1)) +
  ylim(c(-0.2, 0.2)) +
  scale_fill_manual(values = brewer.pal(5, "GnBu")) +
  theme_Publication()

p16.TF <- grid.arrange(p13.TF, p14.TF, p15.TF, ncol = 1)

ggsave(filename = paste0(figure_dir, 'Fig5D.D10.TF.pdf'), p16.TF, width = 3, height = 9)

p13.TF.cl <- p13.TF + 
  xlab('') +
  ylab('') +
  theme(legend.position = "none")

p14.TF.cl <- p14.TF + 
  xlab('') +
  ylab('') +
  theme(legend.position = "none")

p15.TF.cl <- p15.TF +
  xlab('') +
  ylab('') +
  theme(legend.position = "none")

p16.TF.cl <- grid.arrange(p13.TF.cl, p14.TF.cl, p15.TF.cl, ncol = 1)
ggsave(filename = paste0(figure_dir, 'Fig5D.D10.TF.clean.pdf'), p16.TF.cl, width = 2.8, height = 7.5)



##############################################################################################
## Fig 5C
##############################################################################################

# read annotation
dabtram_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/dabtram_saver_cor_SMS.csv')
cocl2_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cocl2_saver_cor_SMS.csv')
cis_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cis_saver_cor_SMS.csv')

# Filter for genes of interest
dabtram_annotations <- dabtram_annotations[dabtram_annotations$Gene.of.interest == 'Yes', ]
cocl2_annotations <- cocl2_annotations[cocl2_annotations$Gene.of.interest == 'Yes', ]
cis_annotations <- cis_annotations[cis_annotations$Gene.of.interest == 'Yes', ]

annotations <- rbind(dabtram_annotations[, c('gene', 'Pathway.')], 
                     cocl2_annotations[, c('gene', 'Pathway.')],
                     cis_annotations[, c('gene', 'Pathway.')]) %>% distinct()
# group annotations by gene and merge pathway information
annotations <- annotations %>% group_by(gene) %>% summarise(Pathway. = paste(Pathway., collapse = ', '))


annotations <- annotations[order(annotations$Pathway.), ]
annotations$order <- seq(1, nrow(annotations))

d0_cor.RNA_anno <- d0_cor.RNA[rownames(d0_cor.RNA) %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]
d10_cor.RNA_anno <- d10_cor.RNA[rownames(d10_cor.RNA) %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]


d0_cor.RNA_anno <- merge(d0_cor.RNA_anno, annotations, by.x = 'row.names', by.y = 'gene')
d0_cor.RNA_anno <- d0_cor.RNA_anno[order(d0_cor.RNA_anno$order), ]
rownames(d0_cor.RNA_anno) <- d0_cor.RNA_anno$Row.names


d10_cor.RNA_anno <- merge(d10_cor.RNA_anno, annotations, by.x = 'row.names', by.y = 'gene')
d10_cor.RNA_anno <- d10_cor.RNA_anno[order(d10_cor.RNA_anno$order), ]
rownames(d10_cor.RNA_anno) <- d10_cor.RNA_anno$Row.names

d0_cor.RNA_anno1 <- d0_cor.RNA_anno[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'Pathway.')]
d0_cor.RNA_anno1$Gene <- rownames(d0_cor.RNA_anno1)

d10_cor.RNA_anno1 <- d10_cor.RNA_anno[, c('correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10', 'Pathway.')]
d10_cor.RNA_anno1$Gene <- rownames(d10_cor.RNA_anno1)

d0_d10.RNA <- merge(d0_cor.RNA_anno1, d10_cor.RNA_anno1, by = c('Pathway.', 'Gene'))
write.csv(d0_d10.RNA, paste0(result_dir, 'd0_d10_correlation_mod.csv'), row.names = F)

d0_d10.RNA <- read.csv(paste0(result_dir, 'd0_d10_correlation_mod.csv'))
rownames(d0_d10.RNA) <- d0_d10.RNA$Gene

pdf(paste0(figure_dir, "Supp_heatmap_genes_d0_d10.pdf"), width = 8, height = 15)
Heatmap(d0_d10.RNA[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = FALSE, row_split = d0_d10.RNA$`Pathway.`, cluster_columns = FALSE,
        show_row_names = T, row_title_rot = 0, 
        column_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()


# read Gavish metaprogram

gavish.mp <- read.csv('/Users/emiliac/Dropbox/Thesis/resources/Signatures/Gavish_Malignant_Meta_Programs.csv')

gavish.mp.list <- list()
for (mp in colnames(gavish.mp)) {
  gavish.mp.list[[mp]] <- c(gavish.mp[[mp]] %>% as.character())
}

gavish.mp.list <- lapply(1:length(gavish.mp.list), function(x) {
  gs <- gavish.mp.list[[x]]
  gs <- unique(na.omit(gs[gs %in% d0_cor.RNA$gene]))
  return(gs)
})
names(gavish.mp.list) <- colnames(gavish.mp)

mps <- c("Cell.Cycle...G2.M", "Cell.Cycle...G1.S", "Cell.Cycle.HMG.rich", "Chromatin",
         "Stress", "Hypoxia", "Stress..in.vitro.", "Proteasomal.degradation", 
         "Unfolded.protein.response", "Protein.maturation", "Translation.initiation", "EMT.I", 
         "EMT.II", "EMT.III", "EMT.IV", "Interferon.MHC.II..I.",
         "Interferon.MHC.II..II.", "Secreted.I", "Secreted.II", "Skin.pigmentation")

gavish.mp.list <- gavish.mp.list[mps]

# turn gavish.mp.list into a dataframe
gavish.mp.df <- data.frame(matrix(ncol = 2, nrow = sum(sapply(gavish.mp.list, length))))
colnames(gavish.mp.df) <- c('gene', 'Gavish_MP')
gavish.mp.df$gene <- unlist(gavish.mp.list)
gavish.mp.df$Gavish_MP <- rep(names(gavish.mp.list), sapply(gavish.mp.list, length))

gavish.mp.df$Gavish_MP_short <- ifelse(grepl('EMT', gavish.mp.df$Gavish_MP), 'EMT', gavish.mp.df$Gavish_MP)
gavish.mp.df$Gavish_MP_short <- ifelse(grepl('Interferon', gavish.mp.df$Gavish_MP_short), 'Interferon', gavish.mp.df$Gavish_MP_short)
gavish.mp.df$Gavish_MP_short <- ifelse(grepl('Cell.Cycle', gavish.mp.df$Gavish_MP_short), 'Cell.Cycle', gavish.mp.df$Gavish_MP_short)
gavish.mp.df$Gavish_MP_short <- ifelse(grepl('Secreted', gavish.mp.df$Gavish_MP_short), 'Secreted', gavish.mp.df$Gavish_MP_short)
gavish.mp.df$Gavish_MP_short <- ifelse(grepl('Stress', gavish.mp.df$Gavish_MP_short), 'Stress', gavish.mp.df$Gavish_MP_short)

d0_d10.RNA <- merge(d0_cor.RNA, d10_cor.RNA, by = 'gene')
d0_d10.RNA <- d0_d10.RNA[, c("gene", "correlation.DABTRAM_d0", "correlation.COCL2_d0", "correlation.CIS_d0", 
                             "correlation.DABTRAM_d10", "correlation.COCL2_d10", "correlation.CIS_d10")]

rownames(d0_d10.RNA) <- d0_d10.RNA$gene

d0_d10.RNA <- merge(d0_d10.RNA, gavish.mp.df, by = 'gene')

d0_d10.RNA$Gavish_MP <- as.factor(gavish.mp.df$Gavish_MP)

pdf(paste0(figure_dir, "Fig5C.heatmap_output_d0_d10_gavish_mp.pdf"), width = 5, height = 7)
Heatmap(as.matrix(d0_d10.RNA[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')]),
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = TRUE, row_split = d0_d10.RNA$`Gavish_MP_short`, cluster_columns = FALSE,
        show_row_names = T, row_title_rot = 0, 
        column_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()

d0_d10.RNA <- t(d0_d10.RNA)
colnames(d0_d10.RNA) <- d0_d10.RNA[1,]
Heatmap(as.matrix(d0_d10.RNA[c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10'), ]),
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = TRUE, column_split = d0_d10.RNA['Gavish_MP_short', ], cluster_columns = FALSE,
        show_row_names = T, row_title_rot = 0, 
        row_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)



## Plot correlation with MP UCEll scores
gavish.mp <- read.csv('/Users/emiliac/Dropbox/Thesis/resources/Signatures/Gavish_Malignant_Meta_Programs.csv')
gavish.mp.names <- as.data.frame(colnames(gavish.mp))
colnames(gavish.mp.names) <- 'Gavish.MP.names'
gavish.mp.names$Gavish.MP.names2 <- paste0('MP', seq(1, nrow(gavish.mp.names)), '-', gavish.mp.names$Gavish.MP.names)

d0_cor.gavish.mp <- read.csv(paste0(result_dir, "GavishMP_UCell_cor_d0_d10.csv"))
d10_cor.gavish.mp <- read.csv(paste0(result_dir, "GavishMP_UCell_cor_d10_w5.csv"))


mps <- c('Cell Cycle - G1/S.rho', 'Cell Cycle - G2/M.rho', 'Cell Cycle HMG-rich.rho', 'Chromatin.rho',
         'Stress.rho', 'Hypoxia.rho', 'Stress (in vitro).rho', 'Proteasomal degradation.rho',
         'Unfolded protein response.rho', 'Protein maturation.rho', 'Translation initiation.rho', 'EMT-I.rho',
         'EMT-II.rho', 'EMT-III.rho', 'EMT-IV.rho', 'Interferon/MHC-II (I).rho',
         'Interferon/MHC-II (II).rho', 'Secreted I.rho', 'Secreted II.rho', 'Skin-pigmentation.rho')

d0_d10_cor.gavish.mp <- merge(d0_cor.gavish.mp, d10_cor.gavish.mp, by = 'MetaProgram')
rownames(d0_d10_cor.gavish.mp) <- d0_d10_cor.gavish.mp$MetaProgram
d0_d10_cor.gavish.mp$MetaProgram <- gsub('.rho', '', d0_d10_cor.gavish.mp$MetaProgram)
d0_d10_cor.gavish.mp$MetaProgram.name <- gsub('\\-', '.', d0_d10_cor.gavish.mp$MetaProgram)
d0_d10_cor.gavish.mp$MetaProgram.name <- gsub(' ', '.', d0_d10_cor.gavish.mp$MetaProgram.name)
d0_d10_cor.gavish.mp$MetaProgram.name <- gsub('\\/', '.', d0_d10_cor.gavish.mp$MetaProgram.name)
d0_d10_cor.gavish.mp$MetaProgram.name <- gsub('\\(', '.', d0_d10_cor.gavish.mp$MetaProgram.name)
d0_d10_cor.gavish.mp$MetaProgram.name <- gsub('\\)', '.', d0_d10_cor.gavish.mp$MetaProgram.name)


d0_d10_cor.gavish.mp <- d0_d10_cor.gavish.mp[mps, ]
d0_d10_cor.gavish.mp <- merge(d0_d10_cor.gavish.mp, gavish.mp.names, by.x = 'MetaProgram.name', by.y = 'Gavish.MP.names', all.x = TRUE)
d0_d10_cor.gavish.mp$Label <- str_split_fixed(d0_d10_cor.gavish.mp$Gavish.MP.names2, '-', 2)[, 1]
d0_d10_cor.gavish.mp$Label <- factor(d0_d10_cor.gavish.mp$Label, levels = paste0('MP', seq(1,41)))
d0_d10_cor.gavish.mp <- d0_d10_cor.gavish.mp[order(d0_d10_cor.gavish.mp$Label), ]

rownames(d0_d10_cor.gavish.mp) <- d0_d10_cor.gavish.mp$Gavish.MP.names2


pdf(paste0(figure_dir, "Fig5C.heatmap_output_d0_d10_gavish_mp.pdf"), width = 5, height = 7)
Heatmap(as.matrix(d0_d10_cor.gavish.mp[, c('DABTRAM.D0.D10', 'COCL2.D0.D10', 'CIS.D0.D10', 'DABTRAM.D10.W5', 'COCL2.D10.W5', 'CIS.D10.W5')]),
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = F, cluster_columns = FALSE,
        show_row_names = T, row_title_rot = 0, 
        column_split = c(rep('D0-to-D10', 3), rep('D10-to-W5', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()

pdf(paste0(figure_dir, "Fig5C.heatmap_output_d0_d10_gavish_mp.pdf"), width = 5, height = 7)
pheatmap(as.matrix(d0_d10_cor.gavish.mp[, c('DABTRAM.D0.D10', 'COCL2.D0.D10', 'CIS.D0.D10', 'DABTRAM.D10.W5', 'COCL2.D10.W5', 'CIS.D10.W5')]),
         col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
         cluster_rows = F,
         cluster_cols = F,
         column_split = c(rep('D0-to-D10', 3), rep('D10-to-W5', 3)),
         cellwidth = 20,
         cellheight = 20,
         border = TRUE)
dev.off()

##############################################################################################
## Fig 5E feature plot
##############################################################################################
remove_unassigned_cells <- TRUE

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day10_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_week5_DABTRAM.RData'))

all_data@misc <- all_data_fatepotential
all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
all_data[["ft.DABTRAM.umap"]] <- all_data_ft_DABTRAM_umap
all_data[["chromVar_day0"]] <- all_data_chromVar_day0
all_data[["chromVar_day10_DABTRAM"]] <- all_data_chromVar_day10_DABTRAM
all_data[["chromVar_week5_DABTRAM"]] <- all_data_chromVar_week5_DABTRAM


# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# ==============================================================================
# wrangle
# ==============================================================================

metadat <- all_data@meta.data
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
metadat.day10_DABTRAM <- metadat[metadat$dataset == 'day10_DABTRAM', ]
metadat.week5_DABTRAM <- metadat[metadat$dataset == 'week5_DABTRAM', ]

# ft umap
ft_umap <- all_data@reductions[["ft.DABTRAM.umap"]]@cell.embeddings
ft_umap <- as.data.frame(ft_umap)
ft_umap$cell_id <- rownames(ft_umap)

# chromVAR
chromVar_day0 <- all_data[["chromVar_day0"]]@data
chromVar_day10_DABTRAM <- all_data[["chromVar_day10_DABTRAM"]]@data
chromVar_week5_DABTRAM <- all_data[["chromVar_week5_DABTRAM"]]@data

chromVar_day0 <- as.data.frame(chromVar_day0[c('JUN'), ])
colnames(chromVar_day0) <- 'JUN'
chromVar_day0$cell_id <- rownames(chromVar_day0)
chromVar_day0 <- chromVar_day0[rownames(metadat.day0), ]

chromVar_day10_DABTRAM <- as.data.frame(chromVar_day10_DABTRAM[c('JUN'), ])
colnames(chromVar_day10_DABTRAM) <- 'JUN'
chromVar_day10_DABTRAM$cell_id <- rownames(chromVar_day10_DABTRAM)
chromVar_day10_DABTRAM <- chromVar_day10_DABTRAM[rownames(metadat.day10_DABTRAM), ]

chromVar_week5_DABTRAM <- as.data.frame(chromVar_week5_DABTRAM[c('JUN'), ])
colnames(chromVar_week5_DABTRAM) <- 'JUN'
chromVar_week5_DABTRAM$cell_id <- rownames(chromVar_week5_DABTRAM)
chromVar_week5_DABTRAM <- chromVar_week5_DABTRAM[rownames(metadat.week5_DABTRAM), ]

chromVar_df <- rbind(chromVar_day0, chromVar_day10_DABTRAM, chromVar_week5_DABTRAM)

ft_umap <- merge(ft_umap, chromVar_df, by = 'cell_id')
ft_umap$JUN <- as.numeric(ft_umap$JUN)
ft_umap$JUN.scale <- scale(ft_umap$JUN)
ft_umap$JUN.scale <- ifelse(ft_umap$JUN.scale < -2, -2, ft_umap$JUN.scale)
ft_umap$JUN.scale <- ifelse(ft_umap$JUN.scale > 2, 2, ft_umap$JUN.scale)

ft_umap$JUN <- ifelse(ft_umap$JUN < -5, -5, ft_umap$JUN)
ft_umap$JUN <- ifelse(ft_umap$JUN > 5, 5, ft_umap$JUN)

p <- ggplot(ft_umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2, color = JUN.scale)) +
  geom_point(size = 0.1) +
  scale_color_gradient2(low = "#604cc3", mid = "bisque", high = "#ffa500", midpoint = 0,  na.value = '#E8E8E8') +
  theme_Publication() +
  theme(legend.position = "right",
        legend.direction = "vertical")
ggsave(paste0(figure_dir, 'SuppFig5E.JUN.fatepotential.pdf'), p, width = 3.5, height = 2)

metadat$cell_id <- rownames(metadat)
ft_umap <- merge(metadat[, c('cell_id', 'fatepotential_DABTRAM_d0_d10', 'fatepotential_DABTRAM_d10_w5')], ft_umap, by = 'cell_id')
p1 <- ggplot(ft_umap, aes(x = fatepotential_DABTRAM_d0_d10, y = `JUN`)) +
  geom_point(size = 0.5) +
  # geom_smooth(method = 'lm', alpha = 0.2) +
  stat_cor(method = 'spearman') +
  ggtitle('DABTRAM (day0)') +
  theme_Publication()

p2 <- ggplot(ft_umap, aes(x = fatepotential_DABTRAM_d10_w5, y = `JUN`)) +
  geom_point(size = 0.5) +
  # geom_smooth(method = 'lm', alpha = 0.2) +
  stat_cor(method = 'spearman', label.x = 0, label.y = -2.5) +
  ggtitle('DABTRAM (day10)') +
  theme_Publication()

p3 <- grid.arrange(p1, p2, ncol = 1)
ggsave(paste0(figure_dir, 'SuppFig5E.JUN.fatepotential.spearman.pdf'), p3, width = 3.5, height = 6.5)
