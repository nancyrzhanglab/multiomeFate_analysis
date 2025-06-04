rm(list = ls())

set.seed(123)

library(Seurat)
library(UCell)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
# out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'

treatment <- 'DABTRAM'

remove_unassigned_cells <- TRUE

date_of_run <- Sys.time()
session_info <- devtools::session_info()

# =============================================================================
# reading data
# =============================================================================

load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))

all_data[['saver']] <- all_data_saver
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

genes <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data_other/Franca_et_al_Resistance_Continuum/Module_Genes.csv',
                  header = TRUE, stringsAsFactors = FALSE)

# turn each column into a list
genes <- lapply(genes, function(x) {
  x <- as.character(x)
  x[x!='']
})

# remove genes that are not in the data
genes <- lapply(genes, function(x) {
  x[x %in% all_data@assays[["saver"]]@var.features]
})


custom.dabtram <- c('ACTB', 'TMEM43', 'TPM4', 'CALM2', 'FN1', 'PALLD', 'LMO7',
                    'ACTN1', 'HSPG2', 'MYOF', 'TNFRSF12A', 'TUBB', 'RCN1', 'CRIM1',
                    'COL5A2', 'SAMD5', 'TPM1', 'OXSR1', 'CBX5')
custom.cis <- c('NUCKS1','TUBB','HMGB2', 'ICMT', 'CBX5', 'TUBA1B', 'ANP32B','TYMS',
                'GMNN', 'USP1', 'NASP', 'TMPO', 'NCAPH', 'TK1', 'TUBG1', 'PRC1',
                'PBK', 'SMC3', 'RRM2', 'RAD51AP1')
custom.cocl2 <- c('GXYLT2', 'ANTXR1', 'CADM1', 'ITGB3', 'BICC1', 'SLC1A4',
                  'CADPS', 'HMGA2', 'TIMP3', 'PTPRG', 'SERPINE2', 'IMMP2L', 'LRMDA',
                  'MFSD12', 'SOX5', 'EPHA3', 'PRKG2', 'IL1RAP', 'SLC44A1', 'KCNQ5')
jackpot = c('SOX10', 'MITF', 'FN1', 'AXL', 'EGFR', 'NT5E',
                    'C1S', 'FRZB', 'SERPINB2', 'SERPINE1', 'NGFR',
          'NDRG1', "FEZF1", 'EGR3', 'VGF',
          'WNT5A', 'POSTN', 'PDGFRB', 'NRG1', 'VEGFC', 'FOSL1',
          'RUNX2', 'LOXL2', 'JUN', 'PDGFRC', 'CD44', 'ID3')
          
DABTRAM = c('AXL', 'EGFP', 'NGFR', 'IGFBP5', 'ANXA1',
                    'IGFBP7', 'JUNB', 'BASP1', 'IER2', 'JUN',
                    'CXCL12', 'ANXA2', 'FOS', 'MMP2', 'GLRX',
                    'IL6ST', 'PRNP', 'FOSB', 'CTSL', 'SLC12A8',
       'TFPI2', 'MYL6', 'IFITM3', 'CAV1', 'CD44')

isg.rs <- c('IFI27', 'IRF7', 'USP18', 'BST2', 'CXCL10',
            'DDX60', 'HERC6', 'HLA-B', 'HLA-G',
            'IFI35', 'IFI44', 'IFI44L', 'IFIT1', 'IFIT3',
            'ISG15', 'LGALS3BP', 'LY6E', 'MX1', 'MX2',
            'OAS3', 'OASL', 'PLSCR1', 'STAT1', 'TRIM14',
            'HSD17B1', 'OAS1', 'CA2', 'CCNA1', 'CXCL1',
            'GALC', 'IFI6', 'IFITM1', 'LAMP3', 'MCL1',
            'ROBO1', 'SLC6A15', 'THBS1', 'TIMP3')

jackpot <- jackpot[jackpot %in% all_data@assays[["saver"]]@var.features]
DABTRAM <- DABTRAM[DABTRAM %in% all_data@assays[["saver"]]@var.features]
isg.rs <- isg.rs[isg.rs %in% all_data@assays[["saver"]]@var.features]



out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
fb.dabtram <- read.csv( paste0(out_dir, 'adapting_bias_thres_0_DABTRAM.csv'))
fb.cis <- read.csv( paste0(out_dir, 'adapting_bias_thres_0_CIS.csv'))
fb.cocl2 <- read.csv( paste0(out_dir, 'adapting_bias_thres_0_COCL2.csv'))

colnames(fb.dabtram)[4] <- 'fate_bias_dabtram'
colnames(fb.cis)[4] <- 'fate_bias_cis'
colnames(fb.cocl2)[4] <- 'fate_bias_cocl2'
# =============================================================================
# Calculate UCell
# =============================================================================
scores <- ScoreSignatures_UCell(all_data@assays[["saver"]]@data, features=genes)
scores.df <- as.data.frame(scores)

scores.df$cell_id <- rownames(scores.df)

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

scores.df <- merge(scores.df, metadat[, c('cell_id', 'dataset')], by = 'cell_id', all.x = TRUE)

scores.df.melt <- melt(scores.df, id.vars = c('cell_id', 'dataset'), variable.name = 'module', value.name = 'score')
scores.df.melt$dataset <- factor(scores.df.melt$dataset, levels = c('day0', 'day10_DABTRAM', 'week5_DABTRAM', 'day10_COCL2', 'week5_COCL2', 'day10_CIS', 'week5_CIS'))
scores.df.melt$module <- factor(scores.df.melt$module, levels = c('Module.A_UCell', 'Module.B_UCell', 'Module.C_UCell', 'Module.D_UCell', 'Module.E_UCell', 'Module.F_UCell'))
ggplot(scores.df.melt, aes(x = dataset, y = score)) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~module, ncol = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = 'Dataset', y = 'UCell score')

ggplot(scores.df.melt, aes(x = module, y = score)) +
  geom_violin(scale = 'width') +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~dataset, ncol = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs( x = 'Module', y = 'UCell score')

# =============================================================================
# Calculate UCell
# =============================================================================
scores <- ScoreSignatures_UCell(all_data@assays[["saver"]]@data, features=list('dabtram.adaptation.genes' = custom.dabtram,
                                                                                 'cis.adaptation.genes' = custom.cis,
                                                                                 'cocl2.adaptation.genes' = custom.cocl2,
                                                                               'jackpot' = jackpot,
                                                                               'dabtram' = DABTRAM,
                                                                               'isg.rs' = isg.rs))
scores.df <- as.data.frame(scores)

scores.df$cell_id <- rownames(scores.df)

metadat <- all_data@meta.data
metadat$cell_id <- rownames(metadat)

scores.df <- merge(scores.df, metadat[, c('cell_id', 'dataset', 
                                          'fatepotential_CIS_d0_d10', 'fatepotential_COCL2_d0_d10', 'fatepotential_DABTRAM_d0_d10')], by = 'cell_id', all.x = TRUE)

scores.df <- merge(scores.df, fb.dabtram[, c('cell_id', 'fate_bias_dabtram')], by = 'cell_id', all.x = TRUE)
scores.df <- merge(scores.df, fb.cis[, c('cell_id', 'fate_bias_cis')], by = 'cell_id', all.x = TRUE)
scores.df <- merge(scores.df, fb.cocl2[, c('cell_id', 'fate_bias_cocl2')], by = 'cell_id', all.x = TRUE)

scores.df.day0 <- scores.df[scores.df$dataset == 'day0', ]

scores.df.day0$fp_dabtram_hi_lo <- ifelse(scores.df.day0$fatepotential_DABTRAM_d0_d10 > quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[4], 'high', 'Other')
scores.df.day0$fp_dabtram_hi_lo <- ifelse(scores.df.day0$fatepotential_DABTRAM_d0_d10 < quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[2], 'low', scores.df.day0$fp_dabtram_hi_lo)
scores.df.day0$fb_dabtram_hi_lo <- ifelse(scores.df.day0$fate_bias_dabtram > quantile(scores.df.day0$fate_bias_dabtram)[4], 'high', 'Other')
scores.df.day0$fb_dabtram_hi_lo <- ifelse(scores.df.day0$fate_bias_dabtram < quantile(scores.df.day0$fate_bias_dabtram)[2], 'low', scores.df.day0$fb_dabtram_hi_lo)
scores.df.day0 <- scores.df.day0[scores.df.day0$fp_dabtram_hi_lo != 'Other' & scores.df.day0$fb_dabtram_hi_lo != 'Other', ]
scores.df.day0$fp_fb <- paste0('FP:', scores.df.day0$fp_dabtram_hi_lo, '_', 'FB:', scores.df.day0$fb_dabtram_hi_lo)

ggplot(scores.df.day0, aes(x = fp_dabtram_hi_lo, y = dabtram.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()

ggplot(scores.df.day0, aes(x = fb_dabtram_hi_lo, y = dabtram.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()

ggplot(scores.df.day0, aes(x = fp_fb, y = dabtram.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  stat_compare_means(comparisons = list(c('FP:high_FB:high', 'FP:high_FB:low'),
                                        c('FP:high_FB:high', 'FP:low_FB:high')), label = 'p.signif') +
  theme_bw()
ggsave(paste0(figure_dir, 'Supp_day0_dabtram_adaptation_genes_UCell_violin.pdf'), 
       width = 4.5, height = 2.5)

ggplot(scores.df.day0, aes(x = fp_fb, y = dabtram.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  stat_compare_means(comparisons = list(c('FP:high_FB:high', 'FP:high_FB:low'),
                                        c('FP:high_FB:high', 'FP:low_FB:high'),
                                        c('FP:high_FB:high', 'FP:low_FB:low')), label = 'p.signif') +
  theme_bw()

scores.df.day0 <- scores.df[scores.df$dataset == 'day0', ]
max_val <- stats::quantile(scores.df.day0$dabtram.adaptation.genes_UCell, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(scores.df.day0$dabtram.adaptation.genes_UCell, 
                           probs = 0.01, 
                           na.rm = TRUE)

scores.df.day0$dabtram.adaptation.genes_UCell_scale <- scales::rescale(
  pmax(pmin(scores.df.day0$dabtram.adaptation.genes_UCell, 
            max_val), 
       min_val),
  to = c(-1, 1)
)
scores.df.day0 <- scores.df.day0[order(scores.df.day0$dabtram.adaptation.genes_UCell_scale), ]
ggplot(scores.df.day0, aes(x = fatepotential_DABTRAM_d0_d10, y = fate_bias_dabtram)) +
  geom_point(aes(color = dabtram.adaptation.genes_UCell_scale)) +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[4], linetype = 'dashed', color = 'black') +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[2], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_dabtram)[4], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_dabtram)[2], linetype = 'dashed', color = 'black') +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', midpoint = 0) +
  theme_bw()
ggsave(paste0(figure_dir, 'Supp_day0_dabtram_adaptation_genes_UCell.pdf'), 
       width = 7, height = 4)

max_val <- stats::quantile(scores.df.day0$jackpot_UCell, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(scores.df.day0$jackpot_UCell, 
                           probs = 0.01, 
                           na.rm = TRUE)

scores.df.day0$jackpot_UCell_scale <- scales::rescale(
  pmax(pmin(scores.df.day0$jackpot_UCell, 
            max_val), 
       min_val),
  to = c(-1, 1)
)
ggplot(scores.df.day0, aes(x = fatepotential_DABTRAM_d0_d10, y = fate_bias_dabtram)) +
  geom_point(aes(color = jackpot_UCell_scale)) +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[4], linetype = 'dashed', color = 'black') +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[2], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_dabtram)[4], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_dabtram)[2], linetype = 'dashed', color = 'black') +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', midpoint = 0) +
  theme_bw()

scores.df.day0 <- scores.df[scores.df$dataset == 'day0', ]
max_val <- stats::quantile(scores.df.day0$isg.rs_UCell, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(scores.df.day0$isg.rs_UCell, 
                           probs = 0.01, 
                           na.rm = TRUE)

scores.df.day0$isg.rs_UCell_scale <- scales::rescale(
  pmax(pmin(scores.df.day0$isg.rs_UCell, 
            max_val), 
       min_val),
  to = c(-1, 1)
)
scores.df.day0 <- scores.df.day0[order(scores.df.day0$isg.rs_UCell_scale), ]
ggplot(scores.df.day0, aes(x = fatepotential_DABTRAM_d0_d10, y = fate_bias_dabtram)) +
  geom_point(aes(color = isg.rs_UCell_scale)) +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[4], linetype = 'dashed', color = 'black') +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_DABTRAM_d0_d10)[2], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_dabtram)[4], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_dabtram)[2], linetype = 'dashed', color = 'black') +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', midpoint = 0) +
  theme_bw()


max_val <- stats::quantile(scores.df.day0$cocl2.adaptation.genes_UCell, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(scores.df.day0$cocl2.adaptation.genes_UCell, 
                           probs = 0.01, 
                           na.rm = TRUE)

scores.df.day0$cocl2.adaptation.genes_UCell_scale <- scales::rescale(
  pmax(pmin(scores.df.day0$cocl2.adaptation.genes_UCell, 
            max_val), 
       min_val),
  to = c(-1, 1)
)
ggplot(scores.df.day0, aes(x = fatepotential_COCL2_d0_d10, y = fate_bias_cocl2)) +
  geom_point(aes(color = cocl2.adaptation.genes_UCell_scale)) +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_COCL2_d0_d10)[4], linetype = 'dashed', color = 'black') +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_COCL2_d0_d10)[2], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_cocl2)[4], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_cocl2)[2], linetype = 'dashed', color = 'black') +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', midpoint = 0) +
  theme_bw()

scores.df.day0 <- scores.df[scores.df$dataset == 'day0', ]
scores.df.day0$fp_cocl2_hi_lo <- ifelse(scores.df.day0$fatepotential_COCL2_d0_d10 > quantile(scores.df.day0$fatepotential_COCL2_d0_d10)[4], 'high', 'Other')
scores.df.day0$fp_cocl2_hi_lo <- ifelse(scores.df.day0$fatepotential_COCL2_d0_d10 < quantile(scores.df.day0$fatepotential_COCL2_d0_d10)[2], 'low', scores.df.day0$fp_cocl2_hi_lo)
scores.df.day0$fb_cocl2_hi_lo <- ifelse(scores.df.day0$fate_bias_cocl2 > quantile(scores.df.day0$fate_bias_cocl2)[4], 'high', 'Other')
scores.df.day0$fb_cocl2_hi_lo <- ifelse(scores.df.day0$fate_bias_cocl2 < quantile(scores.df.day0$fate_bias_cocl2)[2], 'low', scores.df.day0$fb_cocl2_hi_lo)
scores.df.day0 <- scores.df.day0[scores.df.day0$fb_cocl2_hi_lo != 'Other' & scores.df.day0$fp_cocl2_hi_lo != 'Other', ]

# scores.df.day0$fp_cocl2_hi_lo <- ifelse(scores.df.day0$fatepotential_COCL2_d0_d10 > median(scores.df.day0$fatepotential_COCL2_d0_d10), 'high', 'low')
# scores.df.day0$fb_cocl2_hi_lo <- ifelse(scores.df.day0$fate_bias_cocl2 > median(scores.df.day0$fate_bias_cocl2), 'high', 'low')
scores.df.day0$fp_fb <- paste0('FP:', scores.df.day0$fp_cocl2_hi_lo, '_', 'FB:', scores.df.day0$fb_cocl2_hi_lo)

ggplot(scores.df.day0, aes(x = fb_cocl2_hi_lo, y = cocl2.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()

ggplot(scores.df.day0, aes(x = fp_cocl2_hi_lo, y = cocl2.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()
ggplot(scores.df.day0, aes(x = fp_fb, y = cocl2.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()

scores.df.day0 <- scores.df[scores.df$dataset == 'day0', ]
# scores.df.day0$fp_cisplatin_hi_lo <- ifelse(scores.df.day0$fatepotential_CIS_d0_d10 > median(scores.df.day0$fatepotential_CIS_d0_d10), 'high', 'low')
# scores.df.day0$fb_cisplatin_hi_lo <- ifelse(scores.df.day0$fate_bias_cis > median(scores.df.day0$fate_bias_cis), 'high', 'low')
scores.df.day0$fp_cisplatin_hi_lo <- ifelse(scores.df.day0$fatepotential_CIS_d0_d10 > quantile(scores.df.day0$fatepotential_CIS_d0_d10)[4], 'high', 'Other')
scores.df.day0$fp_cisplatin_hi_lo <- ifelse(scores.df.day0$fatepotential_CIS_d0_d10 < quantile(scores.df.day0$fatepotential_CIS_d0_d10)[2], 'low', scores.df.day0$fp_cisplatin_hi_lo)
scores.df.day0$fb_cisplatin_hi_lo <- ifelse(scores.df.day0$fate_bias_cis > quantile(scores.df.day0$fate_bias_cis)[4], 'high', 'Other')
scores.df.day0$fb_cisplatin_hi_lo <- ifelse(scores.df.day0$fate_bias_cis < quantile(scores.df.day0$fate_bias_cis)[2], 'low', scores.df.day0$fb_cisplatin_hi_lo)
scores.df.day0 <- scores.df.day0[scores.df.day0$fb_cisplatin_hi_lo != 'Other' & scores.df.day0$fp_cisplatin_hi_lo != 'Other', ]

scores.df.day0$fp_fb <- paste0('FP:', scores.df.day0$fp_cisplatin_hi_lo, '_', 'FB:', scores.df.day0$fb_cisplatin_hi_lo)

ggplot(scores.df.day0, aes(x = fb_cisplatin_hi_lo, y = cis.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()

ggplot(scores.df.day0, aes(x = fp_cisplatin_hi_lo, y = cis.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()

ggplot(scores.df.day0, aes(x = fp_fb, y = cis.adaptation.genes_UCell)) +
  geom_violin() +
  geom_boxplot() +
  theme_bw()

max_val <- stats::quantile(scores.df.day0$cis.adaptation.genes_UCell, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(scores.df.day0$cis.adaptation.genes_UCell, 
                           probs = 0.01, 
                           na.rm = TRUE)

scores.df.day0$cis.adaptation.genes_UCell_scale <- scales::rescale(
  pmax(pmin(scores.df.day0$cis.adaptation.genes_UCell, 
            max_val), 
       min_val),
  to = c(-1, 1)
)
ggplot(scores.df.day0, aes(x = fatepotential_CIS_d0_d10, y = fate_bias_cis)) +
  geom_point(aes(color = cis.adaptation.genes_UCell_scale)) +
  geom_vline(xintercept = quantile(scores.df.day0$fatepotential_CIS_d0_d10)[3], linetype = 'dashed', color = 'black') +
  # geom_vline(xintercept = quantile(scores.df.day0$fatepotential_CIS_d0_d10)[3], linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = quantile(scores.df.day0$fate_bias_cis)[3], linetype = 'dashed', color = 'black') +
  # geom_hline(yintercept = quantile(scores.df.day0$fate_bias_cis)[2], linetype = 'dashed', color = 'black') +
  scale_color_gradient2(low = 'blue', mid = 'bisque', high = 'red', midpoint = 0) +
  theme_bw()
