rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)

treatment <- 'CIS'
results.dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'

keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
  DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
                   "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
                   "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
                   "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
                   "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")),
  COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
                 "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
                 "CAV1")),
  CIS = sort(c("YY1AP1", "LGALS3", "MCF2L", "TIMM50", "AC207130.1",
               "SLC25A6", "EIF3L", "CTSD", "NQO1", "HNMT", "ZFYVE16",
               "PHACTR1", "TNFRSF14", "RAI14", "TRPM1", "HIST1H1C",
               "HIST2H2AC", "SPARC", "TRIM63", "TUBA1B", "HIST1H1A",
               "HIST1H1D", "PYCARD", "FSTL1", "DCT", "CTSK", "HIST1H4C",
               "GDF15", "HIST1H1B"))
)

# 1) compare gene level output from CoSpar and Fate Bias
# =============================================================================
# Read results
# =============================================================================

# cospar
results.cospar.adapting <- read.csv(paste0(results.dir, 'CoSpar_output/day0_', treatment, '_deg_gene_adapting.csv'), row.names=1)
results.cospar.nonadapting <- read.csv(paste0(results.dir, 'CoSpar_output/day0_', treatment, '_dge_gene_NONadapting.csv'), row.names=1)
results.cospar <- rbind(results.cospar.adapting, results.cospar.nonadapting)

# Fate Bias
load(paste0(results.dir, 'geneSaver_on_day0_cor_vec_', treatment, '.RData'))
results.fatebias <- cor_vec

# =============================================================================
# Wrangle
# =============================================================================
results.cospar <- results.cospar[, -1]
comp_df <- merge(results.cospar, results.fatebias, by='gene')

# =============================================================================
# plot
# =============================================================================
ggplot(comp_df, aes(x=ratio, y=cor)) +
  geom_point() +
  # geom_smooth(method='lm') +
  geom_point(data=comp_df[comp_df$gene %in% keygenes[[treatment]],], color='red', size=3) +
  ggrepel::geom_text_repel(data=comp_df[comp_df$gene %in% keygenes[[treatment]],], 
                           aes(label=gene), color = 'red', box.padding = 0.5) +
  stat_cor() +
  labs(title=paste0('CoSpar vs Fate Bias: ', treatment),
       x='CoSpar: GEX ratio b/w adapting and non-adapting day0 cells',
       y='Cyfer: Correlation b/w GEX and fate bias towards d10-estimated-wk5-adapting') +
  theme_bw()


# 2) compare cell level output from CoSpar and Fate Potential
# =============================================================================
# Read results
# =============================================================================

results.cospar <- read.csv(paste0(results.dir, 'CoSpar_output/day0_', treatment, '_cospar_obs.csv'), row.names=1)
results.fatebias <- read.csv(paste0(results.dir, 'adapting_bias_thres_0_', treatment, '.csv'), row.names=1)

# =============================================================================
# Wrangle
# =============================================================================
results.cospar$cell_id <- rownames(results.cospar)
results.fatebias$cell_id <- rownames(results.fatebias)
comp_df <- merge(results.cospar[, c('cell_id', 'fate_bias_intraclone_transition_map_adapting.non.adapting')], 
                 results.fatebias[, c('cell_id', 'bias')], by='cell_id')

# =============================================================================
# Plot
# =============================================================================

ggplot(comp_df, aes(x=fate_bias_intraclone_transition_map_adapting.non.adapting, y=bias)) +
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor() +
  labs(title=paste0('CoSpar vs Fate Bias: ', treatment),
       x='CoSpar: Intraclone transition map b/w adapting and non-adapting day0 cells',
       y='Fate Bias: Bias towards d10-estimated-wk5-adapting') +
  theme_bw()


