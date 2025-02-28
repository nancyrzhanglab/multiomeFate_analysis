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


dabtram_cor_vec <- read.csv('~/Downloads/cor_vec_TF_DABTRAM.csv')
cocl2_cor_vec <- read.csv('~/Downloads/cor_vec_TF_COCL2.csv')
cis_cor_vec <- read.csv('~/Downloads/cor_vec_TF_CIS.csv')


tfs_all <- dabtram_cor_vec$TF
tfs <- grep('JUN', tfs_all, value = TRUE)
tfs <- c(tfs, grep('FOS', tfs_all, value = TRUE))
tfs <- c(tfs, grep('TEAD', tfs_all, value = TRUE))
tfs <- c(tfs, grep('SNAI', tfs_all, value = TRUE))
tfs <- c(tfs, grep('SOX10', tfs_all, value = TRUE))
tfs <- c(tfs, grep('MITF', tfs_all, value = TRUE))
tfs <- tfs[!grepl('var.2', tfs)]

# tfs <- c(tfs, 'ZEB1', 'TCF4', 'CEBPD', 'CEBPA', 'SOX4')
# ==============================================================================
# Wrangle data
# ==============================================================================

dabtram_cor_vec <- dabtram_cor_vec %>% filter(TF %in% tfs)
cocl2_cor_vec <- cocl2_cor_vec %>% filter(TF %in% tfs)
cis_cor_vec <- cis_cor_vec %>% filter(TF %in% tfs)

dabtram_cor_vec <- dabtram_cor_vec[, c('TF', 'cor')]
cocl2_cor_vec <- cocl2_cor_vec[, c('TF', 'cor')]
cis_cor_vec <- cis_cor_vec[, c('TF', 'cor')]

colnames(dabtram_cor_vec) <- c('TF', 'cor.DABTRAM')
colnames(cocl2_cor_vec) <- c('TF', 'cor.COCL2')
colnames(cis_cor_vec) <- c('TF', 'cor.CIS')

df <- merge(dabtram_cor_vec, cocl2_cor_vec, by = 'TF', all = T)
df <- merge(df, cis_cor_vec, by = 'TF', all = T)

df$TF_fam <- c(rep('AP1', 19), 'MITF', rep('SNAI', 3), 'SOX10', rep('TEAD', 4))
rownames(df) <- df$TF

Heatmap(df[, c('cor.DABTRAM', 'cor.COCL2', 'cor.CIS')],
        cluster_rows = F, show_row_names = T, row_title_rot = 0, row_split = df$TF_fam,
        col = colorRamp2(seq(-0.8,0.8, length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        show_column_names = T, column_title_rot = 0,
        border = TRUE)

