rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE

# keygenes <- list(
#   jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
#                    "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
#                    "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
#                    "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
#                    "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
#   DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
#                    "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
#                    "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
#                    "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
#                    "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44")),
#   COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
#                  "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
#                  "CAV1"))
# )

# keygenes <- list(
#   jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
#                    "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
#                    "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
#                    "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
#                    "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")),
#   DABTRAM = sort(c("AXL", "EGFR", "NGFR", "IGFBP5", "ANXA1",
#                    "IGFBP7", "JUNB", "BASP1", "IER2", "JUN",
#                    "CXCL12", "ANXA2", "FOS", "MMP2", "GLRX",
#                    "IL6ST", "PRNP", "FOSB", "CTSL", "SLC12A8",
#                    "TFPI2", "MYL6", "IFITM3", "CAV1", "CD44"))
# )

keygenes <- list(
  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3"))
)

# keygenes <- list(
#   COCL2 = sort(c("CD44", "FN1", "HPCAL1", "SLC16A3", "IGFBP5",
#                  "COL6A2", "MPC2", "PLIN2", "HLA-A", "IGFBP7",
#                  "CAV1"))
# )
keygenes <- unlist(keygenes)

treatment <- 'CIS'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_saver.RData'))

all_data[['saver']] <- all_data_saver

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# final_fit <- readRDS('~/Downloads/final_fit_d0_w5_DABTRAM_scaled.rds')
# fp_adapting <- as.data.frame(final_fit[["cell_imputed_score"]])

fp_adapting <- readRDS('~/Downloads/cell_imputed_score2_d0_d10_predicted_w5_w_d10_size.CIS_1.rds')
fp_adapting <- as.data.frame(fp_adapting)
colnames(fp_adapting) <- 'FP_adapting'
quantile(fp_adapting$FP_adapting)
quantile(fp_adapting$FP_adapting, probs = seq(0, 1, length.out = 11))

fp_adapting$cell_id <- rownames(fp_adapting)

# =============================================================================
# subset to day0
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')

saver.day0 <- all_data_day0@assays[["saver"]]@scale.data

# =============================================================================
# correlation
# =============================================================================

cor_vec <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(cor_vec) <- c('gene', 'cor', 'p_val')
features <- rownames(saver.day0)
for(g in features) {
  gene.exp <- as.data.frame(saver.day0[g, ])
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
cor_vec$p_val <- as.numeric(cor_vec$p_val)
hist(cor_vec$cor)
cor_vec <- cor_vec %>% arrange(desc(cor))
cor_vec$order <- rownames(cor_vec)
cor_vec$order <- as.numeric(cor_vec$order)
cor_vec$keygene <- ifelse(cor_vec$gene %in% keygenes, 'Key gene', 'Other gene')
cor_vec$p_val_adj <- p.adjust(cor_vec$p_val, method = 'BH')
cor_vec$is_significant <- cor_vec$p_val_adj < 0.05
cor_vec <- cor_vec %>% drop_na()

cor_upper <- cor_vec %>% filter(cor > 0 & is_significant) %>% arrange(desc(cor)) %>% tail(1) %>% pull(cor)
cor_lower <- cor_vec %>% filter(cor < 0 & is_significant) %>% arrange(cor) %>% tail(1) %>% pull(cor)

ggplot(cor_vec, aes(x = order, y = cor)) +
  geom_point() +
  geom_point(data = subset(cor_vec, keygene == 'Key gene'), aes(color = 'key gene')) +
  ggrepel::geom_text_repel(data = subset(cor_vec, keygene == 'Key gene'), aes(label = gene), nudge_y = 0.1) +
  geom_hline(yintercept = cor_upper, linetype = 'dashed', color = 'black') +
  geom_hline(yintercept = cor_lower, linetype = 'dashed', color = 'black') +
  labs(title = paste0('Correlation of gene expression with fate potential (', treatment, ')'),
       x = 'Gene', y = 'Correlation') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


write.csv(cor_vec, paste0('~/Downloads/cor_vec_d0_d10_predicted_w5_w_d10_size_', treatment, '_scaled.csv'), row.names = FALSE)
# =============================================================================
# Compare with day0 cor with d0 to d10 FP
# =============================================================================
# cor_vec <- read.csv(paste0('~/Downloads/cor_vec_', treatment, '_scaled1.csv'))
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
comp_df <- comp_df %>% arrange(keygene)

ggplot(comp_df, aes(x = correlation.d10_w5_FP, y = cor)) +
  geom_point() +
  geom_point(data = subset(comp_df, keygene == 'Key gene'), color = 'red') +
  stat_cor() +
  # geom_smooth(method = 'lm') +
  theme_minimal()
