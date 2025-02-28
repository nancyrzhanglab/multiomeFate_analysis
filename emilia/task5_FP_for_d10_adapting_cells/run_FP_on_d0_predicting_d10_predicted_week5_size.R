rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)
library(ggplot2)

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
keygenes <- unlist(keygenes)

treatment <- 'DABTRAM'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_DABTRAM.RData'))
load(paste0(data_dir, 'Writeup10a_data_peakVI_DABTRAM.RData'))

all_data@misc <- all_data_fatepotential
all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
all_data[['peakVI_DABTRAM']] <- all_data_peakVI_DABTRAM

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

# =============================================================================
# Wrangle data
# =============================================================================

all_data_day10 <- subset(all_data, dataset == paste0('day10_', treatment))
metadat.day10 <- all_data_day10@meta.data
metadat.day10$cell_id <- rownames(metadat.day10)

fp_name <- paste0('fatepotential_', treatment, '_d10_w5')

fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['lineage_imputed_count']])
colnames(fp_d10_w5_treatment) <- c('lineage_imputed_count')
fp_d10_w5_treatment$assigned_lineage <- rownames(fp_d10_w5_treatment)


# =============================================================================
# Prepare input for fate potential
# =============================================================================

# prepare the response variables
all_data_day0 <- subset(all_data, dataset == 'day0')
metadat.day0 <- all_data_day0@meta.data

input_df <- metadat.day0 %>% 
  group_by(assigned_lineage) %>%
  summarize(lineage_current_count = n())

input_df <- left_join(input_df, fp_d10_w5_treatment, by = 'assigned_lineage')
input_df$lineage_imputed_count <- round(input_df$lineage_imputed_count)
input_df$lineage_imputed_count[is.na(input_df$lineage_imputed_count)] <- 0
input_df <- as.data.frame(input_df)
rownames(input_df) <- input_df$assigned_lineage
input_df <- input_df[, -c(1)]

# prepare the input
rna_features <- all_data_day0[["fasttopic_DABTRAM"]]@cell.embeddings
atac_features <- all_data_day0[["peakVI_DABTRAM"]]@cell.embeddings

cell_features <- cbind(rna_features, atac_features)

cell_lineage <- all_data_day0$assigned_lineage
# =============================================================================
# Run fate potential
# =============================================================================
set.seed(123)
colnames(input_df) <- c("now", "future")

lineage_future_count <- input_df$future
names(lineage_future_count) <- rownames(input_df)

colnames(input_df) <- c("day0", "day10_DABTRAM")
tab_mat <- as.matrix(input_df)

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = "day10_DABTRAM",
  lineage_future_count = lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 10,
  tab_mat = tab_mat,
  num_folds = 10,
  verbose = 2
)

final_fit <- multiomeFate:::lineage_cv_finalize(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  fit_res = fit_res,
  lineage_future_count = lineage_future_count
)
lineage_imputed_count <- final_fit$lineage_imputed_count
cell_imputed_score <- final_fit$cell_imputed_score
round(final_fit$coefficient_vec, 2)

hist(cell_imputed_score, breaks = 50)
hist(lineage_imputed_count)

lineage_future_count_df <- as.data.frame(lineage_future_count)
colnames(lineage_future_count_df) <- c('lineage_future_count')
lineage_future_count_df$assigned_lineage <- rownames(lineage_future_count_df)

lineage_imputed_count_df <- as.data.frame(lineage_imputed_count)
colnames(lineage_imputed_count_df) <- c('lineage_imputed_count')
lineage_imputed_count_df$assigned_lineage <- rownames(lineage_imputed_count_df)

comp_df <- merge(lineage_future_count_df, lineage_imputed_count_df, by = 'assigned_lineage')
ggplot(comp_df, aes(x = log10(lineage_future_count+1), y = log10(lineage_imputed_count + 1))) +
  geom_point() +
  stat_cor() +
  # geom_abline(intercept = 0, slope = 1, col = 'red') +
  theme_minimal()

View(cell_imputed_score)

# =============================================================================
# Compare with d0 to d10
# =============================================================================
fp_name2 <- paste0('fatepotential_', treatment, '_d0_d10')
fp_d0_d10_treatment <- as.data.frame(all_data_day10@misc[[fp_name2]][['cell_imputed_score']])

colnames(fp_d0_d10_treatment) <- c('cell_imputed_score_d0_d10')
fp_d0_d10_treatment$cell_id <- rownames(fp_d0_d10_treatment)

cell_imputed_score <- as.data.frame(final_fit$cell_imputed_score)
colnames(cell_imputed_score) <- c('cell_imputed_score_d0.d10_to_week5_FP')
cell_imputed_score$cell_id <- rownames(cell_imputed_score)

comp_df <- merge(fp_d0_d10_treatment, cell_imputed_score, by = 'cell_id')
ggplot(comp_df, aes(x = cell_imputed_score_d0_d10, y = cell_imputed_score_d0.d10_to_week5_FP)) +
  geom_point() +
  stat_cor() +
  theme_minimal()
