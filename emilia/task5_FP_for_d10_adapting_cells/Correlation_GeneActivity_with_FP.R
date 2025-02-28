rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

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
                 "CAV1"))
)

# keygenes <- list(
#   COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
#                  "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
#                  "CAV1"))
# )
keygenes <- unlist(keygenes)

treatment <- 'COCL2'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
all_data_day0 <- subset(all_data, dataset == 'day0')
gene.activities <- readRDS('/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/d0/gene_activities.rds')

final_fit <- readRDS('~/Downloads/final_fit_d0_w5_COCL2_scaled.rds')
fp_adapting <- as.data.frame(final_fit[["cell_imputed_score"]])

# fp_adapting <- readRDS('~/Downloads/cell_imputed_score2_d0_d10_adapting_thres0.5_w_d10_size.CIS_scaled.rds')
fp_adapting <- as.data.frame(fp_adapting)
colnames(fp_adapting) <- 'FP_adapting'
quantile(fp_adapting$FP_adapting)
quantile(fp_adapting$FP_adapting, probs = seq(0, 1, length.out = 11))

fp_adapting$cell_id <- rownames(fp_adapting)

# =============================================================================
# Process gene activities
# =============================================================================
all_data_day0[['gene.activities']] <- CreateAssayObject(counts = gene.activities)
all_data_day0 <- NormalizeData(
  object = all_data_day0,
  assay = 'gene.activities',
  normalization.method = 'LogNormalize',
  scale.factor = median(all_data_day0$nCount_gene.activities)
)
gene.activities.normalized <- all_data_day0@assays[["gene.activities"]]@data

# =============================================================================
# correlation
# =============================================================================

cor_vec <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(cor_vec) <- c('gene', 'cor', 'p_val')
features <- intersect(rownames(d0_saver_cor_vec), rownames(gene.activities))  # rownames(gene.activities)
for(g in features) {
  gene.exp <- as.data.frame(gene.activities[g, ])
  colnames(gene.exp) <- 'gene'
  gene.exp$cell_id <- rownames(gene.exp)
  
  gene.exp <- merge(gene.exp, fp_adapting, by = 'cell_id')
  res <- cor.test(gene.exp$FP_adapting, gene.exp$gene)  
  rho <- res$estimate
  p_val <- res$p.value
  cor_vec <- rbind(cor_vec, c(g, rho, p_val))
}
colnames(cor_vec) <- c('gene', 'cor', 'p_val')
cor_vec$cor <- as.numeric(cor_vec$cor)
hist(cor_vec$cor)
cor_vec <- cor_vec %>% arrange(desc(cor))
cor_vec$order <- rownames(cor_vec)
cor_vec$order <- as.numeric(cor_vec$order)
cor_vec$keygene <- ifelse(cor_vec$gene %in% keygenes, 'Key gene', 'Other gene')

ggplot(cor_vec, aes(x = order, y = cor)) +
  geom_point() +
  geom_point(data = subset(cor_vec, keygene == 'Key gene'), color = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# write.csv(cor_vec, paste0('~/Downloads/cor_vec_', treatment, '_scaled.csv'), row.names = FALSE)
# =============================================================================
# Compare with day0 cor with d0 to d10 FP
# =============================================================================
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
load(paste0(result_dir, 'saver_cor_vec.RData'))
d0_saver_cor_vec <- saver_cor_vec[[paste0(tolower(treatment), '_d0_saver_cor_vec')]]
d0_saver_cor_vec$gene <- rownames(d0_saver_cor_vec)
colnames(d0_saver_cor_vec) <- c('correlation.d0_d10_FP', 'p_value.d0_d10_FP', 'gene')

d10_saver_cor_vec <- saver_cor_vec[[paste0(tolower(treatment), '_d10_saver_cor_vec')]]
d10_saver_cor_vec$gene <- rownames(d10_saver_cor_vec)
colnames(d10_saver_cor_vec) <- c('correlation.d10_w5_FP', 'p_value.d10_w5_FP', 'gene')

comp_df <- merge(d0_saver_cor_vec, cor_vec[, c('gene', 'cor')],  by = 'gene')
comp_df <- merge(d10_saver_cor_vec, comp_df,  by = 'gene')
comp_df$keygene <- ifelse(comp_df$gene %in% keygenes, 'Key gene', 'Other gene')

ggplot(comp_df, aes(x = correlation.d10_w5_FP, y = cor)) +
  geom_point() +
  geom_point(data = subset(comp_df, keygene == 'Key gene'), color = 'red') +
  stat_cor() +
  # geom_smooth(method = 'lm') +
  theme_minimal()
