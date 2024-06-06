library(Seurat)
library(Signac)
library(dplyr)
library(stringr)

# ==============================================================================
# Read data
# ==============================================================================
in_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
day10 <- readRDS(paste0(in_dir, 'Raw_and_Processed/day10_CIS_processed.rds'))
DefaultAssay(day10) <- 'ATAC'
metadat <- day10@meta.data

tp_early <- 'day10'
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

lineage_avg_gp <- metadat %>% 
  group_by(assigned_lineage) %>% 
  summarise(avg_gp = mean(cell_growth_potential),
            clone_size = n())
lineage_avg_gp$isWinner <- ifelse(lineage_avg_gp$avg_gp > 0, 'Winner', 'Non-winner')

metadat <- merge(metadat, lineage_avg_gp, by = 'assigned_lineage')
rownames(metadat) <- metadat$cell_barcode

day10 <- AddMetaData(day10, metadat)

# ==============================================================================
# Get differential peaks
# ==============================================================================
Idents(day10) <- "isWinner"

winner_cells <- metadat[metadat$isWinner == 'Winner', ]$cell_barcode
non_winner_cells <- metadat[metadat$isWinner == 'Non-winner', ]$cell_barcode

da_peaks <- FindMarkers(
  object = day10,
  ident.1 = 'Winner',
  min.pct = 0.2,
  test.use = 'wilcox'
)

da_peaks$neg_log10_p_val <- (-1) * log10(da_peaks$p_val)
ggplot(da_peaks) +
  geom_point(aes(x = avg_log2FC, y = neg_log10_p_val))

# da_peaks_sig <- da_peaks[da_peaks$p_val_adj < 0.05, ]
da_peaks_sig <- da_peaks[da_peaks$p_val < 0.05, ]
da_peaks_sig$region <- rownames(da_peaks_sig)
da_peaks_sig[, c('chr', 'start', 'end')] <- str_split_fixed(da_peaks_sig$region, '-', 3)
da_peaks_sig <- da_peaks_sig[, c('chr', 'start', 'end')]

write.table(da_peaks_sig, '~/Downloads/da_peaks_sig_CIS_no_FDR.bed', sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


