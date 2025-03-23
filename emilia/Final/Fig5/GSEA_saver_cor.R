rm(list = ls())
library(tidyverse)
library(GSEABase)
library(Biobase) #base functions for bioconductor; required by GSEABase
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(enrichplot) # great for making the standard GSEA enrichment plots
library(ggplot2)
library(ggpubr)

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
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="gray90"),
            panel.grid.minor.y = element_line(colour="gray90"),
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


data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig5/'


# ==============================================================================
# Read correlation results
# ==============================================================================
load(paste0(result_dir, 'saver_cor_vec.RData'))

hallmark <- read.gmt(paste0(ref_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome <- read.gmt(paste0(ref_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
threeCA <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))

# ==============================================================================
# GSEA
# ==============================================================================

getAndSortTable <- function(name) {
  cor_vec <- as.data.frame(saver_cor_vec[[name]]) %>% drop_na()
  cor_vec$p.value <- as.numeric(cor_vec$p.value)
  cor_vec$p_adj <- p.adjust(cor_vec$p.value, method = 'BH')
  cor_vec <- cor_vec[cor_vec$p_adj < 0.05, ]
  
  cor_vec$gene <- rownames(cor_vec)
  cor_vec$correlation <- as.numeric(cor_vec$correlation)
  cor_vec <- cor_vec[order(cor_vec$correlation, decreasing = TRUE),]
  
  return(cor_vec)
}

day0_dabtram <- getAndSortTable('dabtram_d0_saver_cor_vec')
day0_cocl2 <- getAndSortTable('cocl2_d0_saver_cor_vec')
day0_cis <- getAndSortTable('cis_d0_saver_cor_vec')

day10_dabtram <- getAndSortTable('dabtram_d10_saver_cor_vec')
day10_cocl2 <- getAndSortTable('cocl2_d10_saver_cor_vec')
day10_cis <- getAndSortTable('cis_d10_saver_cor_vec')

# GSEA DABTRAM
gsea_input.day10.DABTRAM <- day10_dabtram$correlation
names(gsea_input.day10.DABTRAM) <- day10_dabtram$gene

set.seed(123)
GSEA_res.day10.DABTRAM <- GSEA(geneList = gsea_input.day10.DABTRAM, 
                         TERM2GENE = hallmark, 
                         pvalueCutoff = 0.2,
                         seed = T,
                         verbose = F)
GSEA_res.day10.DABTRAM.df <- as_tibble(GSEA_res.day10.DABTRAM@result)
GSEA_res.day10.DABTRAM.df <- GSEA_res.day10.DABTRAM.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.DABTRAM.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM', 'setSize')

# GSEA COCL2
gsea_input.day10.COCL2 <- day10_cocl2$correlation
names(gsea_input.day10.COCL2) <- day10_cocl2$gene

GSEA_res.day10.COCL2 <- GSEA(geneList = gsea_input.day10.COCL2, 
                               TERM2GENE = hallmark, 
                               pvalueCutoff = 0.2,
                               seed = T,
                               verbose = F)
GSEA_res.day10.COCL2.df <- as_tibble(GSEA_res.day10.COCL2@result)
GSEA_res.day10.COCL2.df <- GSEA_res.day10.COCL2.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.COCL2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2', 'setSize')

# GSEA CIS
gsea_input.day10.CIS <- day10_cis$correlation
names(gsea_input.day10.CIS) <- day10_cis$gene

GSEA_res.day10.CIS <- GSEA(geneList = gsea_input.day10.CIS, 
                             TERM2GENE = hallmark, 
                             pvalueCutoff = 0.2,
                             seed = T,
                             verbose = F)
GSEA_res.day10.CIS.df <- as_tibble(GSEA_res.day10.CIS@result)
GSEA_res.day10.CIS.df <- GSEA_res.day10.CIS.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day10.CIS.df) <- c('ID', 'NES.CIS', 'p.adjust.CIS', 'qvalue.CIS', 'setSize')

# ==============================================================================
# Plot
# ==============================================================================

# DABTRAM
GSEA_res.day10.DABTRAM.df <- GSEA_res.day10.DABTRAM.df[order(GSEA_res.day10.DABTRAM.df$NES.DABTRAM, decreasing = TRUE),]
GSEA_res.day10.DABTRAM.df <- GSEA_res.day10.DABTRAM.df[GSEA_res.day10.DABTRAM.df$qvalue.DABTRAM < 0.05,]
GSEA_res.day10.DABTRAM.df$neg_log10_pval.DABTRAM <- -log10(GSEA_res.day10.DABTRAM.df$p.adjust.DABTRAM)

p1 <- ggplot(GSEA_res.day10.DABTRAM.df, aes(x = NES.DABTRAM, y = reorder(ID, NES.DABTRAM), color = neg_log10_pval.DABTRAM, size = setSize)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  xlim(0, 2.5) +
  labs(title = 'GSEA DABTRAM Day 10', x = '', y = '') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# COCL2
GSEA_res.day10.COCL2.df <- GSEA_res.day10.COCL2.df[order(GSEA_res.day10.COCL2.df$NES.COCL2, decreasing = TRUE),]
GSEA_res.day10.COCL2.df <- GSEA_res.day10.COCL2.df[GSEA_res.day10.COCL2.df$qvalue.COCL2 < 0.05,]
GSEA_res.day10.COCL2.df$neg_log10_pval.COCL2 <- -log10(GSEA_res.day10.COCL2.df$p.adjust.COCL2)

p2 <- ggplot(GSEA_res.day10.COCL2.df, aes(x = NES.COCL2, y = reorder(ID, NES.COCL2), color = neg_log10_pval.COCL2, size = setSize)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  xlim(-2, 2.5) +
  labs(title = 'GSEA COCL2 Day 10', x = '', y = '') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# CIS
GSEA_res.day10.CIS.df <- GSEA_res.day10.CIS.df[order(GSEA_res.day10.CIS.df$NES.CIS, decreasing = TRUE),]
GSEA_res.day10.CIS.df <- GSEA_res.day10.CIS.df[GSEA_res.day10.CIS.df$qvalue.CIS < 0.05,]
GSEA_res.day10.CIS.df$neg_log10_pval.CIS <- -log10(GSEA_res.day10.CIS.df$p.adjust.CIS)

p3 <- ggplot(GSEA_res.day10.CIS.df, aes(x = NES.CIS, y = reorder(ID, NES.CIS), color = neg_log10_pval.CIS, size = setSize)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black') +
  xlim(-2, 3) +
  labs(title = 'GSEA CIS Day 10', x = 'Normalized Enrichment Score (NES)', y = '') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p4 <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3, 
                heights = c(16, 9, 10),
          common.legend = T, legend = 'right')

ggsave(paste0(figure_dir, 'SuppFig5.GSEA_day10.pdf'), p4, width = 10, height = 15)

