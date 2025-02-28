rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)

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

# final_fit <- readRDS('~/Downloads/final_fit_d0_w5_COCL2_scaled.rds')
# fp_adapting <- as.data.frame(final_fit[["cell_imputed_score"]])

fp_adapting <- readRDS('~/Downloads/cell_imputed_score2_d0_d10_adapting_w_d10_size.CIS.rds')
fp_adapting <- as.data.frame(fp_adapting)
colnames(fp_adapting) <- 'FP_adapting'
quantile(fp_adapting$FP_adapting)
quantile(fp_adapting$FP_adapting, probs = seq(0, 1, length.out = 11))
fp_adapting$category <- ifelse(fp_adapting$FP_adapting > -0.05956693, 'Adapting_precursor', NA)
fp_adapting$category <- ifelse(fp_adapting$FP_adapting < -0.06321579, 'Non_adapting_precursor', fp_adapting$category)

fp_adapting$cell_id <- rownames(fp_adapting)

table(fp_adapting$category)
# =============================================================================
# subset to day0
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')

saver.day0 <- all_data_day0@assays[["saver"]]@scale.data

adapting_precursors <- fp_adapting %>% 
  filter(category == 'Adapting_precursor') %>%
  select(cell_id) %>% 
  pull()

non_adapting_precursors <- fp_adapting %>% 
  filter(category == 'Non_adapting_precursor') %>%
  select(cell_id) %>% 
  pull()

winner <- as.data.frame(saver.day0['FN1', adapting_precursors])
loser <- as.data.frame(saver.day0['FN1', non_adapting_precursors])

colnames(winner) <- 'FN1'
colnames(loser) <- 'FN1'

winner$category <- 'Adapting_precursor'
loser$category <- 'Non_adapting_precursor'

winner$cell_id <- rownames(winner)
loser$cell_id <- rownames(loser)

comp_df <- rbind(winner, loser)

ggplot(comp_df, aes(x = category, y = FN1)) +
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.2) +
  # geom_jitter(width = 0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cor_vec <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(cor_vec) <- c('gene', 'cor')
features <- rownames(saver.day0)
for(g in features) {
  gene.exp <- as.data.frame(saver.day0[g, ])
  colnames(gene.exp) <- 'gene'
  gene.exp$cell_id <- rownames(gene.exp)
  
  gene.exp <- merge(gene.exp, fp_adapting, by = 'cell_id')
  res <- cor.test(gene.exp$FP_adapting, gene.exp$gene)  
  rho <- res$estimate
  cor_vec <- rbind(cor_vec, c(g, rho))
}
colnames(cor_vec) <- c('gene', 'cor')
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


ggplot(gene.exp, aes(x = FP_adapting, y = CD44)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_minimal()


# =============================================================================
# Compare with day0 cor with d0 to d10 FP
# =============================================================================
result_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task2_correlate_fate_potential_and_features_V2/'
load(paste0(result_dir, 'saver_cor_vec.RData'))
dabtram_d0_saver_cor_vec <- saver_cor_vec[["dabtram_d0_saver_cor_vec"]]
dabtram_d0_saver_cor_vec$gene <- rownames(dabtram_d0_saver_cor_vec)
colnames(dabtram_d0_saver_cor_vec) <- c('correlation.d0_d10_FP', 'p_value.d0_d10_FP', 'gene')

comp_df <- merge(dabtram_d0_saver_cor_vec, cor_vec[, c('gene', 'cor')],  by = 'gene')

ggplot(comp_df, aes(x = correlation.d0_d10_FP, y = cor)) +
  geom_point() +
  stat_cor() +
  geom_smooth(method = 'lm') +
  theme_minimal()

dabtram_d0_saver_cor_vec <- dabtram_d0_saver_cor_vec %>% arrange(desc(correlation.d0_d10_FP))
dabtram_d0_saver_cor_vec$order <- seq(1, nrow(dabtram_d0_saver_cor_vec))
dabtram_d0_saver_cor_vec$keygene <- ifelse(dabtram_d0_saver_cor_vec$gene %in% keygenes, 'Key gene', 'Other gene')

ggplot(dabtram_d0_saver_cor_vec, aes(x = order, y = correlation.d0_d10_FP)) +
  geom_point() +
  geom_point(data = subset(dabtram_d0_saver_cor_vec, keygene == 'Key gene'), color = 'red') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

