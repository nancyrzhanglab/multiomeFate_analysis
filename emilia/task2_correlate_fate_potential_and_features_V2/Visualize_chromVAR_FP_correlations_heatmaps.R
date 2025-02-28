library(tidyverse)
library(ggplot2)
library(GGally)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

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

tfs_all <- rownames(dabtram_d0_chromVAR_cor_vec)
tfs <- grep('JUN', tfs_all, value = TRUE)
tfs <- c(tfs, grep('FOS', tfs_all, value = TRUE))
tfs <- c(tfs, grep('TEAD', tfs_all, value = TRUE))
tfs <- c(tfs, grep('SNAI', tfs_all, value = TRUE))
tfs <- c(tfs, grep('SOX10', tfs_all, value = TRUE))
tfs <- c(tfs, grep('MITF', tfs_all, value = TRUE))
tfs <- tfs[!grepl('var.2', tfs)]

tfs <- c(tfs, 'ZEB1', 'TCF4', 'CEBPD', 'CEBPA', 'SOX4')
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
d0_cor <- merge(dabtram_d0_chromVAR_cor_vec, cocl2_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> dplyr::select(-Row.names)

d0_cor <- merge(d0_cor, cis_d0_chromVAR_cor_vec, by = 'row.names')
rownames(d0_cor) <- d0_cor$Row.names
d0_cor <- d0_cor |> dplyr::select(-Row.names)

# Day10
d10_cor <- merge(dabtram_d10_chromVAR_cor_vec, cocl2_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> dplyr::select(-Row.names)

d10_cor <- merge(d10_cor, cis_d10_chromVAR_cor_vec, by = 'row.names')
rownames(d10_cor) <- d10_cor$Row.names
d10_cor <- d10_cor |> dplyr::select(-Row.names)

d0_cor$TF <- rownames(d0_cor)
d10_cor$TF <- rownames(d10_cor)
d0_d10 <- merge(d0_cor, d10_cor, by = 'TF')
d0_d10 <- d0_d10[d0_d10$TF %in% tfs,]
rownames(d0_d10) <- d0_d10$TF

d0_d10$TF_fam <- c(rep('AP1', 19), 'MITF', rep('SNAI', 3), 'SOX10', rep('TEAD', 4), rep('Others', 5))
pdf("~/Downloads/heatmap_output_d0_d10_TF_all_v2.pdf", width = 5, height = 8)
Heatmap(d0_d10[, c('correlation.DABTRAM_d0', 'correlation.COCL2_d0', 'correlation.CIS_d0', 'correlation.DABTRAM_d10', 'correlation.COCL2_d10', 'correlation.CIS_d10')],
        cluster_rows = T, show_row_names = T, row_title_rot = 0, #row_split = d0_d10$TF_fam,
        cluster_columns = F,
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        column_split = c(rep('Day0', 3), rep('Day10', 3)), show_column_names = T, column_title_rot = 0,
        border = TRUE)
dev.off()
