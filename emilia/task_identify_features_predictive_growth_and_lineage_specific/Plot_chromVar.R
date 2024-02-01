library(Seurat)
library(data.table)
library(ggplot2)
library(ggExtra)

# ==============================================================================
# Read & formatting data
# ==============================================================================
in_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/'
chromVar_dev <- readRDS(paste0(in_dir, 'dev_day10_COCL2_Jun_BS_in_peaks_at_target_vs_nonTarget.rds'))
chromVar_dev_z <- chromVar_dev@assays@data@listData[["z"]]
rownames(chromVar_dev_z) <-  c('Jun_binding_site_in_peaks_targetRNA', 'Jun_binding_site_in_peaks_nonTargetRNA')
chromVar_dev_z <- as.data.frame(t(chromVar_dev_z))
chromVar_dev_z$cell_barcode <- rownames(chromVar_dev_z)
chromVar_dev_z <- melt(chromVar_dev_z, id.vars = 'cell_barcode')

# ==============================================================================
# Plotting
# ==============================================================================
ggplot(chromVar_dev_z, aes(x = variable, y = value)) +
  geom_boxplot() +
  stat_summary(fun=median, geom="point", size=2, color="red")

ggplot(chromVar_dev_z, aes(x = variable, y = value)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", size=2, color="red")

# --------------------------------------------------------------------------------------------------------------

in_dir <- '~/Dropbox/Thesis/Lineage_trace/'

# ==============================================================================
# Read & formatting data
# ==============================================================================
sc <- readRDS(paste0(in_dir, 'data/Shaffer_lab/Raw_and_Processed/day10_DABTRAM_processed.rds'))
chromVar_dev <- readRDS(paste0(in_dir, 'outputs/task4_identify_genes_corr_growth_and_lineage_specific/dev_day10_DABTRAM_Jun_BS_in_peaks_at_target_vs_nonTarget.rds'))
chromVar_dev_z <- chromVar_dev@assays@data@listData[["z"]]
rownames(chromVar_dev_z) <-  c('ap1_peaks_targetRNA', 'ap1_peaks_nonTargetRNA')
chromVar_dev_z <- as.data.frame(t(chromVar_dev_z))
chromVar_dev_z$cell_barcode <- rownames(chromVar_dev_z)

rna_targets <- read.csv(paste0(in_dir, 'outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv'))
rna_targets_non <- read.csv(paste0(in_dir, 'outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv'))

load('~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/day10_DABTRAM/Writeup6n_DABTRAM_day10_lineage-imputation_stepdown_concise-postprocessed.RData')
cell_imputed_count <- cell_imputed_count[!is.na(cell_imputed_count)]
cell_imputed_count <- as.data.frame(cell_imputed_count)
cell_imputed_count$cell_barcode <- rownames(cell_imputed_count)

# ==============================================================================
# Calculate average gene expression from set
# ==============================================================================
DefaultAssay(sc) <- "Saver"
sc <- AddModuleScore(
  object = sc,
  features = list(rna_targets$gene),
  name = 'RNA_targets'
)

sc <- AddModuleScore(
  object = sc,
  features = list(rna_targets_non$gene),
  name = 'non_RNA_targets'
)

rna_set_exp <- sc@meta.data[, c('RNA_targets1', 'non_RNA_targets1')]
rna_set_exp$cell_barcode <- rownames(rna_set_exp)
colnames(rna_set_exp) <- c('target', 'non_target', 'cell_barcode')
rna_set_exp_m <- melt(rna_set_exp, id_vars = 'cell_barcode')
colnames(rna_set_exp_m) <- c("cell_barcode", "variable", "RNA_score" )

colnames(chromVar_dev_z) <- c('target', 'non_target', 'cell_barcode')
chromVar_dev_z_m <- melt(chromVar_dev_z, id_vars = 'cell_barcode')
colnames(chromVar_dev_z_m) <- c("cell_barcode", "variable", "chromVar_score" )

rna_set_exp_m <- merge(rna_set_exp_m, chromVar_dev_z_m, by = c('cell_barcode', 'variable'))
rna_set_exp_m <- merge(rna_set_exp_m, cell_imputed_count, by = 'cell_barcode')

ggplot() +
  geom_point(data = rna_set_exp, aes(x = ap1_peaks_targetRNA, y = RNA_targets1), color = 'red', alpha=0.5) +
  geom_point(data = rna_set_exp, aes(x = ap1_peaks_nonTargetRNA, y = non_RNA_targets1), color = 'blue', alpha=0.5)

p <- ggplot(rna_set_exp_m, aes(x = chromVar_score, y = RNA_score, color=variable)) +
  geom_point(alpha=0.3) +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() + theme(legend.position="bottom")
ggMarginal(p, type="density", alpha=0.5, groupColour = TRUE, groupFill = TRUE)

rna_set_exp_m1 <- rna_set_exp_m[rna_set_exp_m$variable == 'target', ]
rna_set_exp_m2 <- rna_set_exp_m[rna_set_exp_m$variable == 'non_target', ]
ggplot(rna_set_exp_m1, aes(x = chromVar_score, y = RNA_score, color=cell_imputed_count)) +
  geom_point(alpha=0.3) +
  scale_color_gradient2(low = 'blue', mid = 'gray', midpoint = -1, high = 'red', limits = c(-3, 2)) +
  theme_bw() + theme(legend.position="bottom")

p <- ggplot(rna_set_exp_m, aes(x = chromVar_score, y = cell_imputed_count, color = variable)) +
  geom_point(alpha=0.3) + 
  scale_color_manual(values=c("red", "blue")) +
  stat_smooth(geom="line", method = 'lm', alpha=1, lwd = 0.8, se = TRUE) +
  theme_bw() + theme(legend.position="bottom")
ggMarginal(p, type="density", alpha=0.5, groupColour = TRUE, groupFill = TRUE)

cor.test(rna_set_exp_m2$chromVar_score, rna_set_exp_m2$cell_imputed_count)
