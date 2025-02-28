library(tidyverse)
library(ggplot2)
library(GGally)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3"))
)

keygenes <- unlist(keygenes)
# ==============================================================================
# Read data
# ==============================================================================
dabtram_cor_vec <- read.csv('~/Downloads/cor_vec_DABTRAM.csv')
cocl2_cor_vec <- read.csv('~/Downloads/cor_vec_COCL2.csv')
cis_cor_vec <- read.csv('~/Downloads/cor_vec_CIS.csv')

dabtram_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/dabtram_saver_cor_SMS.csv')
cocl2_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cocl2_saver_cor_SMS.csv')
cis_annotations <- read.csv('/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/cis_saver_cor_SMS.csv')

# Filter for genes of interest
dabtram_annotations <- dabtram_annotations[dabtram_annotations$Gene.of.interest == 'Yes', ]
cocl2_annotations <- cocl2_annotations[cocl2_annotations$Gene.of.interest == 'Yes', ]
cis_annotations <- cis_annotations[cis_annotations$Gene.of.interest == 'Yes', ]

annotations <- rbind(dabtram_annotations[, c('gene', 'Pathway.')], 
                     cocl2_annotations[, c('gene', 'Pathway.')],
                     cis_annotations[, c('gene', 'Pathway.')]) %>% distinct()

# group annotations by gene and merge pathway information
annotations <- annotations %>% group_by(gene) %>% summarise(Pathway. = paste(Pathway., collapse = ', '))


annotations <- annotations[order(annotations$Pathway.), ]
annotations$order <- seq(1, nrow(annotations))

# ==============================================================================
# Wrangle data
# ==============================================================================

dabtram_cor_vec <- dabtram_cor_vec[dabtram_cor_vec$gene %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]
cocl2_cor_vec <- cocl2_cor_vec[cocl2_cor_vec$gene %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]
cis_cor_vec <- cis_cor_vec[cis_cor_vec$gene %in% c(dabtram_annotations$gene, cocl2_annotations$gene, cis_annotations$gene), ]

dabtram_cor_vec <- dabtram_cor_vec[, c('gene', 'cor')]
colnames(dabtram_cor_vec) <- c('gene', 'cor.DABTRAM')

cocl2_cor_vec <- cocl2_cor_vec[, c('gene', 'cor')]
colnames(cocl2_cor_vec) <- c('gene', 'cor.COCL2')

cis_cor_vec <- cis_cor_vec[, c('gene', 'cor')]
colnames(cis_cor_vec) <- c('gene', 'cor.CIS')

df <- merge(dabtram_cor_vec, cocl2_cor_vec, by = 'gene')
df <- merge(df, cis_cor_vec, by = 'gene')

# df <- merge(df, annotations, by = 'gene')
df <- df[df$gene %in% keygenes, ]
df <- df[order(df$order), ]
rownames(df) <- df$gene
Heatmap(df[, c('cor.DABTRAM', 'cor.COCL2', 'cor.CIS')],
        col = colorRamp2(seq(-1,1,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)),
        cluster_rows = FALSE,  show_row_names = T, row_title_rot = 0,
        border = TRUE)

