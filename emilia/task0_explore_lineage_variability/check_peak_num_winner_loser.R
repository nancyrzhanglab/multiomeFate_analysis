library(Seurat)
library(Signac)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

# ==============================================================================
# Read data
# ==============================================================================
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
day10 <- readRDS(paste0(in_dir, 'Raw_and_Processed/day0_processed.rds'))
DefaultAssay(day10) <- 'ATAC'
metadat <- day10@meta.data

tp_early <- 'day0'
treatment <- "CIS"
load(paste0(in_dir, "Growth_potential/Writeup6r_", treatment, "_", tp_early, "_lineage-imputation_postprocess.RData"))
dabtram_d10_imputed <- cell_imputed_score[!is.na(cell_imputed_score)]

# ==============================================================================
# Wrangling
# ==============================================================================

dabtram_d10_imputed_df <- data.frame(matrix(nrow=length(dabtram_d10_imputed), ncol = 0))
dabtram_d10_imputed_df$cell_barcode <- names(dabtram_d10_imputed)
dabtram_d10_imputed_df$cell_growth_potential <- dabtram_d10_imputed

metadat$cell_barcode <- rownames(metadat)
metadat <- merge(metadat, dabtram_d10_imputed_df, by = 'cell_barcode')

metadat$isWinner <- ifelse(metadat$cell_growth_potential > 0, 'Winner', 'Non-winner')
rownames(metadat) <- metadat$cell_barcode

# ==============================================================================
# Plottiing
# ==============================================================================
metadat$nFeature_ATAC_norm <- metadat$nFeature_ATAC / metadat$nCount_ATAC
ggplot(metadat, aes(x = nFeature_ATAC_norm, y = cell_growth_potential)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_bw()

cor.test(metadat$nFeature_ATAC_norm, metadat$cell_growth_potential)

# ==================================================================================================
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
day0 <- readRDS(paste0(in_dir, 'Raw_and_Processed/day0_processed.rds'))
day10 <- readRDS(paste0(in_dir, 'Raw_and_Processed/day10_CIS_processed.rds'))
week5 <- readRDS(paste0(in_dir, 'Raw_and_Processed/week5_CIS_processed.rds'))

metadat_day0 <- day0@meta.data
metadat_day10_DABTRAM <- day10@meta.data
metadat_week5_DABTRAM <- week5@meta.data

metadat_day10_COCL2 <- day10@meta.data
metadat_week5_COCL2 <- week5@meta.data

metadat_day10_CIS <- day10@meta.data
metadat_week5_CIS <- week5@meta.data

metadat_day0 <- metadat_day0[, c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC')]

metadat_day10_DABTRAM <- metadat_day10_DABTRAM[, c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC')]
metadat_week5_DABTRAM <- metadat_week5_DABTRAM[, c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC')]

metadat_day10_COCL2 <- metadat_day10_COCL2[, c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC')]
metadat_week5_COCL2 <- metadat_week5_COCL2[, c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC')]

metadat_day10_CIS <- metadat_day10_CIS[, c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC')]
metadat_week5_CIS <- metadat_week5_CIS[, c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC')]

metadat_day0$Sample <- 'Day0'
metadat_day10_DABTRAM$Sample <- 'Day10_DABTRAM'
metadat_week5_DABTRAM$Sample <- 'Week5_DABTRAM'
metadat_day10_COCL2$Sample <- 'Day10_COCL2'
metadat_week5_COCL2$Sample <- 'Week5_COCL2'
metadat_day10_CIS$Sample <- 'Day10_CIS'
metadat_week5_CIS$Sample <- 'Week5_CIS'

to_plot <- rbind(metadat_day0, metadat_day10_DABTRAM)
to_plot <- rbind(to_plot, metadat_week5_DABTRAM)
to_plot <- rbind(to_plot, metadat_day10_COCL2)
to_plot <- rbind(to_plot, metadat_week5_COCL2)
to_plot <- rbind(to_plot, metadat_day10_CIS)
to_plot <- rbind(to_plot, metadat_week5_CIS)
to_plot$nFeature_RNA_norm <- to_plot$nFeature_RNA / to_plot$nCount_RNA
to_plot$nFeature_ATAC_norm <- to_plot$nFeature_ATAC / to_plot$nCount_ATAC

sample_order <- c('Day0', 'Day10_DABTRAM', 'Week5_DABTRAM', 'Day10_COCL2', 'Week5_COCL2', 'Day10_CIS', 'Week5_CIS')
ggplot(to_plot, aes(x = factor(Sample, levels= sample_order), y = nFeature_ATAC)) +
  geom_jitter(width=0.1, alpha=0.2) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5))





