rm(list = ls())

library(clusterProfiler)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(circlize)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '/Users/emiliac/Dropbox/Thesis/resources/GSEA_pathways/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

# ==============================================================================
# Read correlation results
# ==============================================================================
load(paste0(result_dir, 'saver_cor_vec.RData'))

threeCA <- read.gmt(paste0(ref_dir, "c4.3ca.v2024.1.Hs.symbols.gmt"))
# ==============================================================================
# Get and sort table
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

# ==============================================================================
# GSEA
# ==============================================================================

# GSEA DABTRAM
day0_dabtram <- day0_dabtram[!grepl('\\.', day0_dabtram$gene),]
day0_dabtram <- day0_dabtram[order(day0_dabtram$correlation, decreasing = TRUE),]

gsea_input_dabtram <- day0_dabtram$correlation
names(gsea_input_dabtram) <- day0_dabtram$gene

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
day0_cocl2 <- day0_cocl2[!grepl('\\.', day0_cocl2$gene),]
day0_cocl2 <- day0_cocl2[order(day0_cocl2$correlation, decreasing = TRUE),]

gsea_input_cocl2 <- day0_cocl2$correlation
names(gsea_input_cocl2) <- day0_cocl2$gene

GSEA_res_cocl2 <- GSEA(geneList = gsea_input_cocl2, 
                       TERM2GENE = threeCA, 
                       pvalueCutoff = 0.2,
                       seed = T,
                       verbose = F)

GSEA_res_cocl2.df <- as_tibble(GSEA_res_cocl2@result)
GSEA_res_cocl2.df <- GSEA_res_cocl2.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_cocl2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2')

# GSEA CIS
day0_cis <- day0_cis[!grepl('\\.', day0_cis$gene),]
day0_cis <- day0_cis[order(day0_cis$correlation, decreasing = TRUE),]

gsea_input_cis <- day0_cis$correlation
names(gsea_input_cis) <- day0_cis$gene

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
                   show_rownames = T, show_colnames = T, cluster_rows = FALSE, cluster_cols = FALSE, 
                   scale = 'none', border_color = NA, fontsize = 8, cellwidth = 18, cellheight = 18,
                   color = c('#604cc3', '#946983', 'gray', '#de9328', '#ffa500'), # colorRampPalette(brewer.pal(5, "RdPu"))(5)
                   filename = paste0(result_dir, 'GSEA_MPs_day0_heatmap1.pdf'))


write.csv(GSEA_res, paste0(result_dir, 'GSEA_MPs_day0.csv'))



 
