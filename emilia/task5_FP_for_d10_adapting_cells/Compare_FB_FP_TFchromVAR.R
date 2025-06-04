rm(list = ls())

set.seed(123)

library(Seurat)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
# out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

treatment <- 'DABTRAM'

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_chromVar_day0.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data[['chromVar.day0']] <- all_data_chromVar_day0
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

out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
fb.dabtram <- read.csv( paste0(out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
fb.cis <- read.csv( paste0(out_dir, 'adapting_bias_thres_0_CIS.csv'))
fb.cocl2 <- read.csv( paste0(out_dir, 'adapting_bias_thres_0_COCL2.csv'))

colnames(fb.dabtram)[4] <- 'fate_bias_dabtram'
colnames(fb.cis)[4] <- 'fate_bias_cis'
colnames(fb.cocl2)[4] <- 'fate_bias_cocl2'

# =============================================================================
# Wrangle data
# =============================================================================
metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

chromVAR_cur <- all_data@assays[["chromVar.day0"]]@data
chromVAR_cur <- as.matrix(chromVAR_cur)

chromVAR_cur.subset <- chromVAR_cur[, c('JUNB', 'FOS', 'TEAD1', 'CTCF', 'JUN', 'HIF1A')]
chromVAR_cur.subset <- as.data.frame(chromVAR_cur.subset)
chromVAR_cur.subset$cell_id <- rownames(chromVAR_cur.subset)

chromVAR_cur.subset <- merge(chromVAR_cur.subset, metadat[, c('cell_id', 'dataset', 
                                          'fatepotential_CIS_d0_d10', 'fatepotential_COCL2_d0_d10', 'fatepotential_DABTRAM_d0_d10')], by = 'cell_id', all.x = TRUE)
chromVAR_cur.subset <- merge(chromVAR_cur.subset, fb.dabtram[, c('cell_id', 'fate_bias_dabtram')], by = 'cell_id', all.x = TRUE)

chromVAR_cur.subset <- chromVAR_cur.subset[chromVAR_cur.subset$dataset == 'day0', ]

chromVAR_cur.subset$fp_dabtram_hi_lo <- ifelse(chromVAR_cur.subset$fatepotential_DABTRAM_d0_d10 > quantile(chromVAR_cur.subset$fatepotential_DABTRAM_d0_d10)[4], 'high', 'Other')
chromVAR_cur.subset$fp_dabtram_hi_lo <- ifelse(chromVAR_cur.subset$fatepotential_DABTRAM_d0_d10 < quantile(chromVAR_cur.subset$fatepotential_DABTRAM_d0_d10)[2], 'low', chromVAR_cur.subset$fp_dabtram_hi_lo)
chromVAR_cur.subset$fb_dabtram_hi_lo <- ifelse(chromVAR_cur.subset$fate_bias_dabtram > quantile(chromVAR_cur.subset$fate_bias_dabtram)[4], 'high', 'Other')
chromVAR_cur.subset$fb_dabtram_hi_lo <- ifelse(chromVAR_cur.subset$fate_bias_dabtram < quantile(chromVAR_cur.subset$fate_bias_dabtram)[2], 'low', chromVAR_cur.subset$fb_dabtram_hi_lo)
chromVAR_cur.subset <- chromVAR_cur.subset[chromVAR_cur.subset$fp_dabtram_hi_lo != 'Other' & chromVAR_cur.subset$fb_dabtram_hi_lo != 'Other', ]
chromVAR_cur.subset$fp_fb <- paste0('FP:', chromVAR_cur.subset$fp_dabtram_hi_lo, '_', 'FB:', chromVAR_cur.subset$fb_dabtram_hi_lo)

chromVAR_cur.subset.melt <- reshape2::melt(chromVAR_cur.subset[, c('cell_id', 'fp_fb', 'JUNB', 'FOS', 'TEAD1', 'CTCF', 'JUN', 'HIF1A')], id.vars = c('cell_id', 'fp_fb'))

ggplot(chromVAR_cur.subset, aes(x = fp_fb, y = HIF1A)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  stat_compare_means(comparisons = list(c('FP:high_FB:high', 'FP:high_FB:low'),
                                        c('FP:high_FB:high', 'FP:low_FB:high')), label = 'p.signif') +
  theme_bw()

ggplot(chromVAR_cur.subset.melt, aes(x = fp_fb, y = value)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5) +
  stat_compare_means(comparisons = list(c('FP:high_FB:high', 'FP:high_FB:low'),
                                        c('FP:high_FB:high', 'FP:low_FB:high')), label = 'p.signif') +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

