rm(list = ls())
library(tidyverse)
library(data.table)
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
d0_cor <- merge(dabtram_d0_saver_cor_vec, cocl2_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

d0_cor <- merge(d0_cor, cis_d0_saver_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> select(-Row.names)

d0_cor$gene <- rownames(d0_cor)

# make a master dataframe for Day10
d10_cor <- merge(dabtram_d10_saver_cor_vec, cocl2_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

d10_cor <- merge(d10_cor, cis_d10_saver_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> select(-Row.names)

d10_cor$gene <- rownames(d10_cor)

# =============================================================================
# Plot density
# =============================================================================

# Plot params
x_min <- -1
x_max <- 1


# Define key genes
res_genes <- c("WNT5A", "AXL", "EGFR", "PDGFRB", "JUN", "NGFR", "PCNA")

# Pick out res genes
keygenes.dabtram.d0 <- d0_cor[d0_cor$gene %in% res_genes, 'correlation.DABTRAM_d0']
keygenes.cocl2.d0 <- d0_cor[d0_cor$gene %in% res_genes, 'correlation.COCL2_d0']
keygenes.cis.d0 <- d0_cor[d0_cor$gene %in% res_genes, 'correlation.CIS_d0']

keygenes.dabtram.d10 <- d10_cor[d10_cor$gene %in% res_genes, 'correlation.DABTRAM_d10']
keygenes.cocl2.d10 <- d10_cor[d10_cor$gene %in% res_genes, 'correlation.COCL2_d10']
keygenes.cis.d10 <- d10_cor[d10_cor$gene %in% res_genes, 'correlation.CIS_d10']

# DABTRAM day0
d0_cor <- d0_cor %>% arrange(correlation.DABTRAM_d0)
d0_cor$order.DABTRAM <- 1:nrow(d0_cor)
p1 <- ggplot(d0_cor, aes(x = order.DABTRAM, y = correlation.DABTRAM_d0)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = d0_cor[d0_cor$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, gene %in% res_genes), 
                           aes(label = gene)) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(title = 'day0 correlation', ylab = 'Correlation', xlab = 'Gene rank')

# COCL2 day0
d0_cor <- d0_cor %>% arrange(correlation.COCL2_d0)
d0_cor$order.COCL2 <- 1:nrow(d0_cor)
p2 <- ggplot(d0_cor, aes(x = order.COCL2, y = correlation.COCL2_d0)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = d0_cor[d0_cor$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, gene %in% res_genes), 
                           aes(label = gene)) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

# CIS day0
d0_cor <- d0_cor %>% arrange(correlation.CIS_d0)
d0_cor$order.CIS <- 1:nrow(d0_cor)
p3 <- ggplot(d0_cor, aes(x = order.CIS, y = correlation.CIS_d0)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = d0_cor[d0_cor$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d0_cor, gene %in% res_genes), 
                           aes(label = gene)) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

p4 <- grid.arrange(p1, p2, p3, ncol = 1)




# DABTRAM day10
d10_cor <- d10_cor %>% arrange(correlation.DABTRAM_d10)
d10_cor$order.DABTRAM <- 1:nrow(d10_cor)
p5 <- ggplot(d10_cor, aes(x = order.DABTRAM, y = correlation.DABTRAM_d10)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = d10_cor[d10_cor$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, gene %in% res_genes), 
                           aes(label = gene)) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

# COCL2 day10
d10_cor <- d10_cor %>% arrange(correlation.COCL2_d10)
d10_cor$order.COCL2 <- 1:nrow(d10_cor)
p6 <- ggplot(d10_cor, aes(x = order.COCL2, y = correlation.COCL2_d10)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = d10_cor[d10_cor$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, gene %in% res_genes), 
                           aes(label = gene)) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')


# CIS day10
d10_cor <- d10_cor %>% arrange(correlation.CIS_d10)
d10_cor$order.CIS <- 1:nrow(d10_cor)
p7 <- ggplot(d10_cor, aes(x = order.CIS, y = correlation.CIS_d10)) +
  geom_point(size = 1, color = 'gray') +
  geom_point(data = d10_cor[d10_cor$gene %in% res_genes, ], color = 'red') +
  ggrepel::geom_text_repel(data = subset(d10_cor, gene %in% res_genes), 
                           aes(label = gene)) +
  ylim(-1, 1) +
  theme_Publication() +
  labs(ylab = 'Correlation', xlab = 'Gene rank')

p8 <- grid.arrange(p5, p6, p7, ncol = 1)

ggsave(filename = paste0(figure_dir, 'Fig5D.pdf'), p8, width = 4, height = 5)

##############################################################################################
## Fig 5D: Comparison of correlations
##############################################################################################

# Day0 correlations
p9 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p10 <- ggplot(d0_cor, aes(x = correlation.CIS_d0, y = correlation.COCL2_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()


p11 <- ggplot(d0_cor, aes(x = correlation.DABTRAM_d0, y = correlation.CIS_d0)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.5, label.y = 0.45) +
  xlim(c(-0.505, 0.505)) +
  ylim(c(-0.505, 0.505)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p12 <- grid.arrange(p9, p10, p11, ncol = 3)

# D10 correlations
p13 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p14 <- ggplot(d10_cor, aes(x = correlation.CIS_d10, y = correlation.COCL2_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()


p15 <- ggplot(d10_cor, aes(x = correlation.DABTRAM_d10, y = correlation.CIS_d10)) +
  # geom_point(color = 'pink') +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor(method = "spearman", label.x = -0.95, label.y = 0.8) +
  xlim(c(-1, 1)) +
  ylim(c(-1, 1)) +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

p16 <- grid.arrange(p13, p14, p15, ncol = 3)

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

fatepotential.all_data <- metadat.all_data[, c('cell_id', 'dataset',
                                               'fatepotential_CIS_d0_d10', 'fatepotential_CIS_d10_w5', 
                                               'fatepotential_COCL2_d0_d10', 'fatepotential_COCL2_d10_w5',
                                               'fatepotential_DABTRAM_d0_d10', 'fatepotential_DABTRAM_d10_w5')]

cor.to.plot <- merge(fatepotential.all_data, saver.all_data.toplot, by = 'cell_id')

# =============================================================================
# Plot data
# =============================================================================
cor.to.plot.melt <- melt(cor.to.plot, id.vars = c('cell_id', 'dataset', 'fatepotential_CIS_d0_d10', 'fatepotential_CIS_d10_w5', 
                                                  'fatepotential_COCL2_d0_d10', 'fatepotential_COCL2_d10_w5',
                                                  'fatepotential_DABTRAM_d0_d10', 'fatepotential_DABTRAM_d10_w5'))
colnames(cor.to.plot.melt)[9] <- 'Gene'
colnames(cor.to.plot.melt)[10] <- 'Saver'
cor.to.plot.melt <- melt(cor.to.plot.melt, id.vars = c('cell_id', 'dataset', 'Gene', 'Saver'))
colnames(cor.to.plot.melt)[5] <- 'FatePotential.index'
colnames(cor.to.plot.melt)[6] <- 'FatePotential.estimate'

cor.to.plot.melt <- cor.to.plot.melt %>% drop_na()
ggplot(cor.to.plot.melt, aes(x = Saver, y = FatePotential.estimate)) +
  geom_point(color = 'gray', size = 1) +
  stat_cor(method = "spearman", label.sep = "\n", na.rm = T, label.x = -0.95, label.y = 0) +
  facet_grid(FatePotential.index ~ Gene, scales = 'free') +
  theme_bw()



