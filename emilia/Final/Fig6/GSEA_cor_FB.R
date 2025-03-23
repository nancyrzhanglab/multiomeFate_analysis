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
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

hallmark <- read.gmt(paste0(ref_dir, "h.all.v2024.1.Hs.symbols.gmt"))
reactome <- read.gmt(paste0(ref_dir, "c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
threeCA <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))
# ==============================================================================
# Read correlation results
# ==============================================================================
load(paste0(result_dir, 'geneSaver_on_day0_cor_vec_DABTRAM.RData'))
cor_vec.DABTRAM <- cor_vec

load(paste0(result_dir, 'geneSaver_on_day0_cor_vec_COCL2.RData'))
cor_vec.COCL2 <- cor_vec

load(paste0(result_dir, 'geneSaver_on_day0_cor_vec_CIS.RData'))
cor_vec.CIS <- cor_vec

# ==============================================================================
# GSEA
# ==============================================================================

getAndSortTable <- function(cor_vec) {
  cor_vec$p_val <- as.numeric(cor_vec$p_val)
  cor_vec$p_adj <- p.adjust(cor_vec$p_val, method = 'BH')
  cor_vec <- cor_vec[cor_vec$p_adj < 0.05, ]
  
  cor_vec$cor <- as.numeric(cor_vec$cor)
  cor_vec <- cor_vec[order(cor_vec$cor, decreasing = TRUE),]
  
  return(cor_vec)
}

cor_vec.DABTRAM <- getAndSortTable(cor_vec.DABTRAM)
cor_vec.COCL2 <- getAndSortTable(cor_vec.COCL2)
cor_vec.CIS <- getAndSortTable(cor_vec.CIS)


# GSEA DABTRAM
gsea_input.day0.DABTRAM <- cor_vec.DABTRAM$cor
names(gsea_input.day0.DABTRAM) <- cor_vec.DABTRAM$gene

set.seed(123)
GSEA_res.day0.DABTRAM <- GSEA(geneList = gsea_input.day0.DABTRAM, 
                               TERM2GENE = hallmark, 
                               pvalueCutoff = 0.2,
                               seed = T,
                               verbose = F)
GSEA_res.day0.DABTRAM.df <- as_tibble(GSEA_res.day0.DABTRAM@result)
GSEA_res.day0.DABTRAM.df <- GSEA_res.day0.DABTRAM.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day0.DABTRAM.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM', 'setSize')

# GSEA COCL2
gsea_input.day0.COCL2 <- cor_vec.COCL2$cor
names(gsea_input.day0.COCL2) <- cor_vec.COCL2$gene

GSEA_res.day0.COCL2 <- GSEA(geneList = gsea_input.day0.COCL2, 
                             TERM2GENE = hallmark, 
                             pvalueCutoff = 0.2,
                             seed = T,
                             verbose = F)
GSEA_res.day0.COCL2.df <- as_tibble(GSEA_res.day0.COCL2@result)
GSEA_res.day0.COCL2.df <- GSEA_res.day0.COCL2.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day0.COCL2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2', 'setSize')

# GSEA CIS
gsea_input.day0.CIS <- cor_vec.CIS$cor
names(gsea_input.day0.CIS) <- cor_vec.CIS$gene

GSEA_res.day0.CIS <- GSEA(geneList = gsea_input.day0.CIS, 
                           TERM2GENE = hallmark, 
                           pvalueCutoff = 0.2,
                           seed = T,
                           verbose = F)
GSEA_res.day0.CIS.df <- as_tibble(GSEA_res.day0.CIS@result)
GSEA_res.day0.CIS.df <- GSEA_res.day0.CIS.df[, c('ID', 'NES', 'p.adjust', 'qvalue', 'setSize')]
colnames(GSEA_res.day0.CIS.df) <- c('ID', 'NES.CIS', 'p.adjust.CIS', 'qvalue.CIS', 'setSize')

# merge
GSEA_res <- merge(GSEA_res.day0.DABTRAM.df, GSEA_res.day0.COCL2.df, by = c('ID'), all = T)
GSEA_res <- merge(GSEA_res, GSEA_res.day0.CIS.df, by = c('ID'), all = T)
GSEA_res <- GSEA_res[, c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'NES.COCL2', 'p.adjust.COCL2', 'NES.CIS', 'p.adjust.CIS')]
GSEA_res_m <- melt(GSEA_res, id.vars = 'ID')

# filter out ones with two or mor NAs
GSEA_res <- GSEA_res[rowSums(is.na(GSEA_res)) <= 4,]
GSEA_res.DABTRAM <- GSEA_res[, c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM')]
GSEA_res.COCL2 <- GSEA_res[, c('ID', 'NES.COCL2', 'p.adjust.COCL2')]
GSEA_res.CIS <- GSEA_res[, c('ID', 'NES.CIS', 'p.adjust.CIS')]

colnames(GSEA_res.DABTRAM) <- c('ID', 'NES', 'p.adjust')
colnames(GSEA_res.COCL2) <- c('ID', 'NES', 'p.adjust')
colnames(GSEA_res.CIS) <- c('ID', 'NES', 'p.adjust')

GSEA_res.DABTRAM$treatment <- 'DABTRAM'
GSEA_res.COCL2$treatment <- 'COCL2'
GSEA_res.CIS$treatment <- 'CIS'

GSEA_res.to_plot <- rbind(GSEA_res.DABTRAM, GSEA_res.COCL2, GSEA_res.CIS)
GSEA_res.to_plot$p.adjust <- ifelse(GSEA_res.to_plot$p.adjust < 0.05, GSEA_res.to_plot$p.adjust, NA)
GSEA_res.to_plot$neg_log10_pval <- -log10(GSEA_res.to_plot$p.adjust)
GSEA_res.to_plot$NES.abs <- abs(GSEA_res.to_plot$NES)
GSEA_res.to_plot$treatment <- factor(GSEA_res.to_plot$treatment, levels = c('DABTRAM', 'COCL2', 'CIS'))

# plot bubble plot
ggplot(GSEA_res.to_plot, aes(x = treatment, y = ID, size = NES.abs, color = NES)) +
  geom_point() +
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  labs(title = 'GSEA', x = '', y = '') +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(paste0(figure_dir, 'GSEA_day0_FB_cor.pdf'), width = 7, height = 8)
       