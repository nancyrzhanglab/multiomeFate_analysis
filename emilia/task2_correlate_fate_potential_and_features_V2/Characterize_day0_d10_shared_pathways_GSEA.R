library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(circlize)

result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

GSEA_res.d0 <- read.csv(paste0(result_dir, 'GSEA_MPs_day0.csv'), row.names = 1)
GSEA_res.d10 <- read.csv(paste0(result_dir, 'GSEA_MPs_day10.csv'), row.names = 1)

GSEA_res.d0 <- GSEA_res.d0[ , c('ID', 'NES.DABTRAM', 'NES.COCL2', 'NES.CIS', 'ID_short')]
GSEA_res.d10 <- GSEA_res.d10[ , c('ID', 'NES.DABTRAM', 'NES.COCL2', 'NES.CIS', 'ID_short')]

colnames(GSEA_res.d0) <- c('ID', 'NES.DABTRAM.d0', 'NES.COCL2.d0', 'NES.CIS.d0', 'ID_short')
colnames(GSEA_res.d10) <- c('ID', 'NES.DABTRAM.d10', 'NES.COCL2.d10', 'NES.CIS.d10', 'ID_short')

GSEA_res <- merge(GSEA_res.d0, GSEA_res.d10, by = c('ID', 'ID_short'), all = T)
rownames(GSEA_res) <- GSEA_res$ID_short
order <- c('MP2_CELL_CYCLE_G1_S', 'MP4_CHROMATIN', 'MP32_SKIN_PIGMENTATION', 
           'MP5_STRESS', 'MP6_HYPOXIA', 'MP7_STRESS_IN_VITRO', 'MP12_EMT_1', 'MP13_EMT_2', 'MP14_EMT_3', 'MP15_EMT_4', 'MP39_METAL_RESPONSE',
           'MP17_INTERFERON_MHC_II_1', 'MP18_INTERFERON_MHC_II_2', 'MP22_SECRETED_1', 'MP3_CELL_CYLCE_HMG_RICH', 
           'MP8_PROTEASOMAL_DEGRADATION', 'MP9_UNFOLDED_PROTEIN_RESPONSE', 'MP10_PROTEIN_MATURATION')
GSEA_res <- GSEA_res[order, ]
pheatmap::pheatmap(as.matrix(GSEA_res[, c('NES.DABTRAM.d0', 'NES.COCL2.d0', 'NES.CIS.d0',
                                          'NES.DABTRAM.d10', 'NES.COCL2.d10', 'NES.CIS.d10')]), 
                   show_rownames = T, show_colnames = T, cluster_rows = FALSE, cluster_cols = FALSE, 
                   scale = 'none', border_color = 'black', fontsize = 8, cellwidth = 18, cellheight = 18,
                   gaps_col = c(3, 6), 
                   color = c('#604cc3', '#946983', 'gray', '#de9328', '#ffa500'))
