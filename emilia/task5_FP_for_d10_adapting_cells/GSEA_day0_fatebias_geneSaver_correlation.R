rm(list = ls())

library(clusterProfiler)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(circlize)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

# ==============================================================================
# Read correlation results
# ==============================================================================
load(paste0(out_dir, 'geneSaver_on_day0_cor_vec_DABTRAM.RData'))
cor_vec.DABTRAM <- cor_vec

load(paste0(out_dir, 'geneSaver_on_day0_cor_vec_COCL2.RData'))
cor_vec.COCL2 <- cor_vec

load(paste0(out_dir, 'geneSaver_on_day0_cor_vec_CIS.RData'))
cor_vec.CIS <- cor_vec


threeCA <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))
# ==============================================================================
# Get and sort table
# ==============================================================================

getAndSortTable <- function(df) {
  cor_vec <- as.data.frame(df) %>% drop_na()
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

# ==============================================================================
# GSEA
# ==============================================================================

# GSEA DABTRAM
cor_vec.DABTRAM <- cor_vec.DABTRAM[!grepl('\\.', cor_vec.DABTRAM$gene),]
cor_vec.DABTRAM <- cor_vec.DABTRAM[order(cor_vec.DABTRAM$cor, decreasing = TRUE),]

gsea_input_dabtram <- cor_vec.DABTRAM$cor
names(gsea_input_dabtram) <- cor_vec.DABTRAM$gene

set.seed(123)
GSEA_res_dabtram <- GSEA(geneList = gsea_input_dabtram, 
                         TERM2GENE = threeCA, 
                         pvalueCutoff = 0.2,
                         seed = T,
                         verbose = F)
GSEA_res_dabtram.df <- as_tibble(GSEA_res_dabtram@result)
GSEA_res_dabtram.df <- GSEA_res_dabtram.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_dabtram.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM')

# GSEA COCL2
cor_vec.COCL2 <- cor_vec.COCL2[!grepl('\\.', cor_vec.COCL2$gene),]
cor_vec.COCL2 <- cor_vec.COCL2[order(cor_vec.COCL2$cor, decreasing = TRUE),]

gsea_input_cocl2 <- cor_vec.COCL2$cor
names(gsea_input_cocl2) <- cor_vec.COCL2$gene

GSEA_res_cocl2 <- GSEA(geneList = gsea_input_cocl2, 
                       TERM2GENE = threeCA, 
                       pvalueCutoff = 0.2,
                       seed = T,
                       verbose = F)

GSEA_res_cocl2.df <- as_tibble(GSEA_res_cocl2@result)
GSEA_res_cocl2.df <- GSEA_res_cocl2.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_cocl2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2')

# GSEA CIS
cor_vec.CIS <- cor_vec.CIS[!grepl('\\.', cor_vec.CIS$gene),]
cor_vec.CIS <- cor_vec.CIS[order(cor_vec.CIS$cor, decreasing = TRUE),]

gsea_input_cis <- cor_vec.CIS$cor
names(gsea_input_cis) <- cor_vec.CIS$gene

GSEA_res_cis <- GSEA(geneList = gsea_input_cis, 
                     TERM2GENE = threeCA, 
                     pvalueCutoff = 0.2,
                     seed = T,
                     verbose = F)

GSEA_res_cis.df <- as_tibble(GSEA_res_cis@result)
GSEA_res_cis.df <- GSEA_res_cis.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_cis.df) <- c('ID', 'NES.CIS', 'p.adjust.CIS', 'qvalue.CIS')

# plot

GSEA_res <- merge(GSEA_res_dabtram.df, GSEA_res_cocl2.df, by = 'ID', all = T)
GSEA_res <- merge(GSEA_res, GSEA_res_cis.df, by = 'ID', all = T)
GSEA_res$n_NA <- rowSums(is.na(GSEA_res[, c('NES.DABTRAM', 'NES.COCL2', 'NES.CIS')]))
GSEA_res <- GSEA_res[GSEA_res$n_NA < 2,]

GSEA_res <- GSEA_res[order(GSEA_res$n_NA, GSEA_res$NES.DABTRAM, decreasing = c(FALSE, TRUE)),]
GSEA_res$ID_short <- gsub('GAVISH_3CA_MALIGNANT_METAPROGRAM_', 'MP', GSEA_res$ID)
rownames(GSEA_res) <- GSEA_res$ID_short
pheatmap::pheatmap(as.matrix(GSEA_res[, c('NES.DABTRAM', 'NES.COCL2', 'NES.CIS')]), 
                   show_rownames = T, show_colnames = T, cluster_rows = TRUE, cluster_cols = FALSE, 
                   scale = 'none', border_color ='black', fontsize = 8, cellwidth = 18, cellheight = 18,
                   color = c('#604cc3', '#946983', 'gray', '#de9328', '#ffa500'), # colorRampPalette(brewer.pal(5, "RdPu"))(5)
                   filename = paste0(out_dir, 'GSEA_MPs_day0_heatmap.pdf'))
