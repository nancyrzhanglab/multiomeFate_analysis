library(Seurat)
library(data.table)
library(ggplot2)
library(ggExtra)

# ==============================================================================
# Read & formatting data
# ==============================================================================
in_dir <- '~/Dropbox/Thesis/Lineage_trace/'

sc <- readRDS(paste0(in_dir, 'data/Shaffer_lab/Raw_and_Processed/day10_DABTRAM_processed.rds'))
rna_targets <- read.csv(paste0(in_dir, 'outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2.csv'))
rna_targets_non <- read.csv(paste0(in_dir, 'outputs/task4_identify_genes_corr_growth_and_lineage_specific/common_genes_in_corr_with_growth_v2_nonTarget.csv'))

dev <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/ChromVar/Vierstra/mat.motifs_day10_DABTRAM_Vierstra_annotations.rds")
dev_z <- as.data.frame(dev@assays@data@listData[["z"]])

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


motifs <- rownames(dev_z)[grep("FOS|JUN|TEAD", rownames(dev_z))]
dev_z <- as.data.frame(t(dev_z))
dev_z <- dev_z[, motifs]
dev_z$cell_barcode <- rownames(dev_z)

rna_tf <- merge(rna_set_exp_m, dev_z, by = 'cell_barcode')
rna_tf <- merge(rna_tf, cell_imputed_count, by = 'cell_barcode')
# ==============================================================================
# Calculate average gene expression from set
# ==============================================================================
p <- ggplot(rna_tf, aes(x = `AC0201|TEAD|TEA`, y = RNA_score, col = variable)) +
  geom_point(alpha = 0.3) +
  stat_smooth(geom="line", method = 'lm', alpha=1, lwd = 0.8) +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() + theme(legend.position="bottom")
ggMarginal(p, type="density", alpha=0.5, groupColour = TRUE, groupFill = TRUE)


rna_tf_sub <- rna_tf[rna_tf$variable == 'target', ]
ggplot(rna_tf_sub, aes(x = `AC0242|FOSL/JUND|bZIP`, y = RNA_score, col = cell_imputed_count)) +
  geom_point(alpha = 0.3) +
  scale_color_gradient2(low = 'blue', mid = 'gray', midpoint = -1, high = 'red', limits = c(-3, 2)) +
  theme_bw()
  
ggplot(rna_tf_sub, aes(x = cell_imputed_count, y = RNA_score)) +
  geom_point(alpha = 0.3) +
  scale_color_gradient2(low = 'blue', mid = 'gray', midpoint = -1, high = 'red', limits = c(-3, 2)) +
  theme_bw()

