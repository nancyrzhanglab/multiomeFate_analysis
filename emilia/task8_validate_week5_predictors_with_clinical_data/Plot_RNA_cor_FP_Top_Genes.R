rm(list = rm())
library(Seurat)
library(UCell)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(GGally)
library(ggdensity)
library(RColorBrewer)
library(hdrcde)
library(ggdensity)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
ref_dir <- '~/Downloads/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'

# Custom signatures
custom.cis <- c('NUCKS1','TUBB','HMGB2', 'ICMT', 'CBX5', 'TUBA1B', 'ANP32B','TYMS',
                'GMNN', 'USP1', 'NASP', 'TMPO', 'NCAPH', 'TK1', 'TUBG1', 'PRC1',
                'PBK', 'SMC3', 'RRM2', 'RAD51AP1')
custom.dabtram <- c('ACTB', 'TMEM43', 'TPM4', 'CALM2', 'FN1', 'PALLD', 'LMO7',
                    'ACTN1', 'HSPG2', 'MYOF', 'TNFRSF12A', 'TUBB', 'RCN1', 'CRIM1',
                    'COL5A2', 'SAMD5', 'TPM1', 'OXSR1', 'CBX5')
custom.cocl2 <- c('GXYLT2', 'ANTXR1', 'CADM1', 'ITGB3', 'BICC1', 'SLC1A4',
                  'CADPS', 'HMGA2', 'TIMP3', 'PTPRG', 'SERPINE2', 'IMMP2L', 'LRMDA',
                  'MFSD12', 'SOX5', 'EPHA3', 'PRKG2', 'IL1RAP', 'SLC44A1', 'KCNQ5')

# ==============================================================================
# Read data general
# ==============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data@misc <- all_data_fatepotential
all_data[["Saver"]] <- all_data_saver

genes.all <- all_data_saver@data@Dimnames[[1]]

custom.list <- list(custom.dabtram, custom.cocl2, custom.cis)
custom.list <- lapply(1:length(custom.list), function(x) {
  gs <- custom.list[[x]]
  gs <- unique(na.omit(gs[gs %in% genes.all]))
  return(gs)
})
names(custom.list) <- c('custom.dabtram', 'custom.cocl2', 'custom.cis')

# ==============================================================================
# Calcualte module scores 
# ==============================================================================

metadat <- all_data@meta.data
scores <- ScoreSignatures_UCell(all_data@assays[["Saver"]]@data, 
                                features=custom.list)

colnames(scores) <- names(custom.list)
scores.df <- as.data.frame(scores) 

# ==============================================================================
# Compare with fate potentials
# ==============================================================================
fatepot.DABTRAM.d0.d10 <- all_data_fatepotential[['fatepotential_DABTRAM_d0_d10']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.DABTRAM.d0.d10) <- 'fatepotential_DABTRAM_d0_d10'

fatepot.DABTRAM.d10.w5 <- all_data_fatepotential[['fatepotential_DABTRAM_d10_w5']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.DABTRAM.d10.w5) <- 'fatepotential_DABTRAM_d10_w5'
fatepot.DABTRAM.d10.w5$cell_id <- rownames(fatepot.DABTRAM.d10.w5)

fatepot.COCL2.d0.d10 <- all_data_fatepotential[['fatepotential_COCL2_d0_d10']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.COCL2.d0.d10) <- 'fatepotential_COCL2_d0_d10'
fatepot.COCL2.d0.d10$cell_id <- rownames(fatepot.COCL2.d0.d10)

fatepot.COCL2.d10.w5 <- all_data_fatepotential[['fatepotential_COCL2_d10_w5']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.COCL2.d10.w5) <- 'fatepotential_COCL2_d10_w5'
fatepot.COCL2.d10.w5$cell_id <- rownames(fatepot.COCL2.d10.w5)

fatepot.CIS.d0.d10 <- all_data_fatepotential[['fatepotential_CIS_d0_d10']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.CIS.d0.d10) <- 'fatepotential_CIS_d0_d10'
fatepot.CIS.d0.d10$cell_id <- rownames(fatepot.CIS.d0.d10)

fatepot.CIS.d10.w5 <- all_data_fatepotential[['fatepotential_CIS_d10_w5']][["cell_imputed_score"]] %>% as.data.frame()
colnames(fatepot.CIS.d10.w5) <- 'fatepotential_CIS_d10_w5'
fatepot.CIS.d10.w5$cell_id <- rownames(fatepot.CIS.d10.w5)

# add fate potential
scores.df <- merge(scores.df, fatepot.DABTRAM.d0.d10, by='row.names', all = T)
colnames(scores.df)[1] <- 'cell_id'
scores.df <- merge(scores.df, fatepot.DABTRAM.d10.w5, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.COCL2.d0.d10, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.COCL2.d10.w5, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.CIS.d0.d10, by='cell_id', all = T)
scores.df <- merge(scores.df, fatepot.CIS.d10.w5, by='cell_id', all = T)

# ==============================================================================
# Subset to day10
# ==============================================================================
metadat.day10.cis <- metadat[metadat$dataset == 'day10_CIS', ]
scores.df.day10.cis <- scores.df[scores.df$cell_id %in% rownames(metadat.day10.cis), ]

fatepot.d10.w5 <- c('fatepotential_DABTRAM_d10_w5', 'fatepotential_COCL2_d10_w5', 'fatepotential_CIS_d0_d10')
scores.df.day10.cis <- scores.df.day10.cis[, -which(names(scores.df.day10.cis) %in% fatepot.d10.w5)]

# ==============================================================================
# Subset to day0
# ==============================================================================
metadat.day0 <- metadat[metadat$dataset == 'day0', ]
scores.df.day0 <- scores.df[scores.df$cell_id %in% rownames(metadat.day0), ]

fatepot.d0.d10 <- c('fatepotential_DABTRAM_d10_w5', 'fatepotential_COCL2_d10_w5', 'fatepotential_CIS_d10_w5')
scores.df.day0 <- scores.df.day0[, -which(names(scores.df.day0) %in% fatepot.d0.d10)]

# ==============================================================================
# Plot
# ==============================================================================

ggplot(scores.df.day10.cis, aes(x=fatepotential_CIS_d10_w5, y=custom.cis)) +
  # geom_point() +
  geom_hdr(aes(fill = after_stat(probs)), color = "black", alpha = 0.8, probs = c(0.99, 0.8, 0.6, 0.4, 0.2)) +
  stat_cor() +
  labs( y='CIS w5 top20 predictors') +
  scale_fill_manual(values = brewer.pal(5, "RdPu")) +
  theme_Publication()

ggplot(scores.df.day0, aes(x=fatepotential_CIS_d0_d10, y=custom.cis)) +
  geom_point() +
  stat_cor(method = 'spearman') +
  labs( y='CIS w5 top20 predictors') +
  theme_Publication()
