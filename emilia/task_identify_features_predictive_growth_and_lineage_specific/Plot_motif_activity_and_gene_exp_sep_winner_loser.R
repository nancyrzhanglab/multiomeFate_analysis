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
day10 <- readRDS(paste0(in_dir, 'Raw_and_Processed/day10_COCL2_processed.rds'))
DefaultAssay(day10) <- 'ATAC'
metadat <- day10@meta.data

tp_early <- 'day10'
treatment <- "COCL2"
load(paste0(in_dir, "Growth_potential/Writeup6r_", treatment, "_", tp_early, "_lineage-imputation_postprocess.RData"))
dabtram_d10_imputed <- cell_imputed_score[!is.na(cell_imputed_score)]

load(paste0(in_dir, '/ChromVar/JASPAR/chromVar_day10_data.RData'))
motif_names <- read.csv(paste0(in_dir, 'motif_info.csv'))

rna_targets <- read.csv('~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv')
rna_targets_non <- read.csv('~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv')

# ==============================================================================
# Wrangling
# ==============================================================================

dabtram_d10_imputed_df <- data.frame(matrix(nrow=length(dabtram_d10_imputed), ncol = 0))
dabtram_d10_imputed_df$cell_barcode <- names(dabtram_d10_imputed)
dabtram_d10_imputed_df$cell_growth_potential <- dabtram_d10_imputed

metadat$cell_barcode <- rownames(metadat)
metadat <- merge(metadat, dabtram_d10_imputed_df, by = 'cell_barcode')

# lineage_avg_gp <- metadat %>% 
#   group_by(assigned_lineage) %>% 
#   summarise(avg_gp = mean(cell_growth_potential),
#             clone_size = n())
# lineage_avg_gp$isWinner <- ifelse(lineage_avg_gp$avg_gp > 0, 'Winner', 'Non-winner')

# metadat <- merge(metadat, lineage_avg_gp, by = 'assigned_lineage')
metadat$isWinner <- ifelse(metadat$cell_growth_potential > 0, 'Winner', 'Non-winner')
rownames(metadat) <- metadat$cell_barcode

day10 <- AddMetaData(day10, metadat)

chromvar_results_dabtram <- as.data.frame(chromvar_results_cocl2)
chromvar_results_dabtram$motif_code <- rownames(chromvar_results_dabtram)
chromvar_results_dabtram <- merge(chromvar_results_dabtram, motif_names, by = 'motif_code')

chromvar_results_dabtram_use <- chromvar_results_dabtram[chromvar_results_dabtram$motif_names %in% c('JUN', 'FOS::JUN', 
                                                                                                     "TEAD1", "RORA", "RXRG"), ]
# ==============================================================================
# Calculate average gene expression from set
# ==============================================================================
DefaultAssay(day10) <- "Saver"
day10 <- AddModuleScore(
  object = day10,
  features = list(rna_targets$gene),
  name = 'RNA_targets'
)

day10 <- AddModuleScore(
  object = day10,
  features = list(rna_targets_non$gene),
  name = 'non_RNA_targets'
)

rna_set_exp <- day10@meta.data[, c('RNA_targets1', 'non_RNA_targets1')]
rna_set_exp$cell_barcode <- rownames(rna_set_exp)
colnames(rna_set_exp) <- c('target', 'non_target', 'cell_barcode')
rna_set_exp_m <- melt(rna_set_exp, id_vars = 'cell_barcode')
colnames(rna_set_exp_m) <- c("cell_barcode", "variable", "RNA_score" )

chromvar_results_dabtram_use <- t(chromvar_results_dabtram_use)
chromvar_results_dabtram_use <- chromvar_results_dabtram_use[-c(1), ]
colnames(chromvar_results_dabtram_use) <- chromvar_results_dabtram_use[nrow(chromvar_results_dabtram_use), ]
chromvar_results_dabtram_use <- chromvar_results_dabtram_use[-c(nrow(chromvar_results_dabtram_use)), ]
chromvar_results_dabtram_use <- as.data.frame(chromvar_results_dabtram_use)
chromvar_results_dabtram_use$cell_barcode <- rownames(chromvar_results_dabtram_use)

rna_set_exp <- merge(rna_set_exp, chromvar_results_dabtram_use, by = 'cell_barcode')
rna_set_exp <- merge(rna_set_exp, metadat[, c('cell_barcode', 'isWinner')], by = 'cell_barcode')

rna_set_exp$RORA <- as.numeric(rna_set_exp$RORA)
rna_set_exp$RXRG <- as.numeric(rna_set_exp$RXRG)
rna_set_exp$TEAD1 <- as.numeric(rna_set_exp$TEAD1)
rna_set_exp$`FOS::JUN` <- as.numeric(rna_set_exp$`FOS::JUN`)
rna_set_exp$`JUN` <- as.numeric(rna_set_exp$`JUN`)
ggplot(rna_set_exp) +
  geom_point(aes(x = `JUN`, y = target)) +
  facet_wrap(. ~ isWinner)

ggplot(rna_set_exp, aes(x = isWinner, y = non_target)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  theme_bw()

ggplot(rna_set_exp, aes(x = isWinner, y = TEAD1)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  theme_bw()
ggplot(rna_set_exp, aes(x = `JUN`, y = target, col = isWinner)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = 'lm') +
  theme_bw()


rna_set_exp_m <- pivot_longer(rna_set_exp, cols = c("RORA", "TEAD1", "FOS::JUN", "JUN", "RXRG"), 
                              names_to = 'TF', names_prefix = '', values_to = 'ChromVAR')

ggplot(rna_set_exp_m, aes(x = ChromVAR, y = target)) +
  geom_point(alpha=0.5, aes(color= isWinner)) +
  scale_color_manual(values = c('#0079FF', '#FF004D')) +
  geom_smooth(method = 'lm', aes(linetype = isWinner), color = 'black') +
  facet_wrap(. ~ TF, scale='free') +
  theme_bw()

ggplot(rna_set_exp_m, aes(x = ChromVAR, y = non_target, col = isWinner)) +
  geom_point(alpha=0.5, aes(color= isWinner)) +
  scale_color_manual(values = c('#0079FF', '#FF004D')) +
  geom_smooth(method = 'lm', aes(linetype = isWinner), color = 'black') +
  facet_wrap(. ~ TF, scale='free') +
  theme_bw()

