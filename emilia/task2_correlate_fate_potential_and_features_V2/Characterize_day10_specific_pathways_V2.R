library(Seurat)
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


day10_dabtram <- getAndSortTable('dabtram_d10_saver_cor_vec')
day10_cocl2 <- getAndSortTable('cocl2_d10_saver_cor_vec')
day10_cis <- getAndSortTable('cis_d10_saver_cor_vec')


# ==============================================================================
# GSEA
# ==============================================================================

# GSEA DABTRAM
day10_dabtram <- day10_dabtram[!grepl('\\.', day10_dabtram$gene),]
day10_dabtram <- day10_dabtram[order(day10_dabtram$correlation, decreasing = TRUE),]

gsea_input_dabtram <- day10_dabtram$correlation
names(gsea_input_dabtram) <- day10_dabtram$gene

set.seed(123)
GSEA_res_dabtram <- GSEA(geneList = gsea_input_dabtram, 
                         TERM2GENE = threeCA, 
                         pvalueCutoff = 0.1,
                         seed = T,
                         verbose = F)
GSEA_res_dabtram.df <- as_tibble(GSEA_res_dabtram@result)
GSEA_res_dabtram.df <- GSEA_res_dabtram.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_dabtram.df) <- c('ID', 'NES.DABTRAM', 'p.adjust.DABTRAM', 'qvalue.DABTRAM')

# GSEA COCL2
day10_cocl2 <- day10_cocl2[!grepl('\\.', day10_cocl2$gene),]
day10_cocl2 <- day10_cocl2[order(day10_cocl2$correlation, decreasing = TRUE),]

gsea_input_cocl2 <- day10_cocl2$correlation
names(gsea_input_cocl2) <- day10_cocl2$gene

GSEA_res_cocl2 <- GSEA(geneList = gsea_input_cocl2, 
                       TERM2GENE = threeCA, 
                       pvalueCutoff = 0.2,
                       seed = T,
                       verbose = F)

GSEA_res_cocl2.df <- as_tibble(GSEA_res_cocl2@result)
GSEA_res_cocl2.df <- GSEA_res_cocl2.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_cocl2.df) <- c('ID', 'NES.COCL2', 'p.adjust.COCL2', 'qvalue.COCL2')

# GSEA CIS
day10_cis <- day10_cis[!grepl('\\.', day10_cis$gene),]
day10_cis <- day10_cis[order(day10_cis$correlation, decreasing = TRUE),]

gsea_input_cis <- day10_cis$correlation
names(gsea_input_cis) <- day10_cis$gene

GSEA_res_cis <- GSEA(geneList = gsea_input_cis, 
                     TERM2GENE = threeCA, 
                     pvalueCutoff = 0.2,
                     seed = T,
                     verbose = F)

GSEA_res_cis.df <- as_tibble(GSEA_res_cis@result)
GSEA_res_cis.df <- GSEA_res_cis.df[, c('ID', 'NES', 'p.adjust', 'qvalue')]
colnames(GSEA_res_cis.df) <- c('ID', 'NES.CIS', 'p.adjust.CIS', 'qvalue.CIS')

# ==============================================================================
# Merge GSEA results
# ==============================================================================
GSEA_res <- merge(GSEA_res_dabtram.df, GSEA_res_cocl2.df, by = 'ID', all = T)
GSEA_res <- merge(GSEA_res, GSEA_res_cis.df, by = 'ID', all = T)

GSEA_res$n_NA <- rowSums(is.na(GSEA_res[, c('NES.DABTRAM', 'NES.COCL2', 'NES.CIS')]))
# GSEA_res <- GSEA_res[GSEA_res$n_NA < 2,]

GSEA_res <- GSEA_res[order(GSEA_res$NES.DABTRAM, GSEA_res$NES.COCL2, GSEA_res$NES.CIS, decreasing = T),]
GSEA_res$ID_short <- gsub('GAVISH_3CA_MALIGNANT_METAPROGRAM_', 'MP', GSEA_res$ID)

rownames(GSEA_res) <- GSEA_res$ID_short
pheatmap::pheatmap(as.matrix(GSEA_res[, c('NES.DABTRAM', 'NES.COCL2', 'NES.CIS')]), 
                   show_rownames = T, show_colnames = T, cluster_rows = FALSE, cluster_cols = FALSE, 
                   scale = 'none', border_color = NA, fontsize = 8, cellwidth = 18, cellheight = 18,
                   color = c('#604cc3', '#946983', 'gray', '#de9328', '#ffa500'),
                   # color = c("#604cc3", "#654fbd", "#6a51b8", "#6e54b2", "#7357ac", 
                   #           "#785aa6", "#7c5ca1", "#805e9c", "#846097", "#886292", 
                   #           "#8c658d", "#906788", "#936984", "#976b80", "#9a6d7b", 
                   #           "#9e6e77", "#a17073", "#a4726f", "#a8746b", "#ab7668", 
                   #           "#ae7864", "#b17960", "#b47b5c", "#b87d58", "#bb7e54", 'gray',
                   #           "#be8050", "#c0824d", "#c38349", "#c68546", "#c98742", 
                   #           "#cc883f", "#cf8a3c", "#d18c38", "#d48d35", "#d78f32", 
                   #           "#da902e", "#dc922b", "#df9328", "#e29524", "#e49621", 
                   #           "#e7981e", "#ea991a", "#ec9b17", "#ef9c14", "#f29e10", 
                   #           "#f49f0d", "#f7a10a", "#faa207", "#fca403", "#ffa500"),
                   filename = paste0(result_dir, 'GSEA_MPs_day10_heatmap2.pdf'))

write.csv(GSEA_res, paste0(result_dir, 'GSEA_MPs_day10.csv'))

