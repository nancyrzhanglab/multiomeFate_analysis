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
keygenes <- unlist(keygenes)

treatment <- 'COCL2'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_', treatment, '.RData'))
load(paste0(data_dir, 'Writeup10a_data_peakVI_', treatment, '.RData'))

all_data@misc <- all_data_fatepotential
# all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
# all_data[['peakVI_DABTRAM']] <- all_data_peakVI_DABTRAM

all_data[['fasttopic_COCL2']] <- all_data_fasttopic_COCL2
all_data[['peakVI_COCL2']] <- all_data_peakVI_COCL2

# all_data[['fasttopic_CIS']] <- all_data_fasttopic_CIS
# all_data[['peakVI_CIS']] <- all_data_peakVI_CIS


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

all_data_week5 <- subset(all_data, dataset == paste0('week5_', treatment))
metadat.week5 <- all_data_week5@meta.data
metadat.week5$cell_id <- rownames(metadat.week5)

lineage_size <- metadat.week5 %>% 
  group_by(assigned_lineage) %>%
  summarize(n_cells = n())

quantile(lineage_size$n_cells)
hist(lineage_size$n_cells, breaks = 50)


# =============================================================================
# Prepare input for fate potential
# =============================================================================

# prepare the response variables
all_data_day0 <- subset(all_data, dataset == 'day0')
metadat.day0 <- all_data_day0@meta.data

input_df <- metadat.day0 %>% 
  group_by(assigned_lineage) %>%
  summarize(lineage_current_count = n())


input_df <- merge(input_df, lineage_size, by = 'assigned_lineage', all.x = TRUE)
input_df[is.na(input_df)] <- 0


ggplot(input_df, aes(x = lineage_current_count, y = n_cells)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = 'red') +
  theme_minimal()


cell_lineage <- all_data_day0$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))

# prepare the input
rna_features <- all_data_day0[[paste0("fasttopic_", treatment)]]@cell.embeddings
atac_features <- all_data_day0[[paste0("peakVI_", treatment)]]@cell.embeddings

cell_features <- cbind(rna_features, atac_features)
# cell_features <- scale(cell_features)

# =============================================================================
# Run fate potential
# =============================================================================
set.seed(123)

rownames(input_df) <- input_df$assigned_lineage
input_df <- input_df[, c('lineage_current_count', 'n_cells')]
input_df <- as.matrix(input_df)
colnames(input_df) <- c("now", "future")

lineage_future_count <- input_df[, 'future']

colnames(input_df) <- c("day0", paste0("week5_", treatment))
tab_mat <- input_df

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = paste0("week5_", treatment),
  lineage_future_count = lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 10,
  tab_mat = tab_mat,
  num_folds = 20,
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
max <- log10(max(comp_df$lineage_future_count, comp_df$lineage_imputed_count)) + 0.13
ggplot(comp_df, aes(x = log10(lineage_future_count+1), y = log10(lineage_imputed_count + 1))) +
  geom_jitter(width = 0.1) +
  stat_cor() +
  geom_abline(intercept = 0, slope = 1, col = 'red') +
  # lims(x = c(0, max), y = c(0, max)) +
  labs(x = 'log10(lineage_week5_count + 1)', y = 'log10(day0-predicted-week5-lineage-size + 1)', title = treatment) +
  theme_bw()

View(cell_imputed_score)

saveRDS(final_fit, '~/Downloads/final_fit_d0_w5_COCL2_1.rds')


# =============================================================================
# Compare with d0 to d10
# =============================================================================
fp_name2 <- paste0('fatepotential_', treatment, '_d0_d10')
fp_d0_d10_treatment <- as.data.frame(all_data_day10@misc[[fp_name2]][['cell_imputed_score']])

colnames(fp_d0_d10_treatment) <- c('cell_imputed_score_d0_d10')
fp_d0_d10_treatment$cell_id <- rownames(fp_d0_d10_treatment)

cell_imputed_score <- as.data.frame(final_fit$cell_imputed_score)
colnames(cell_imputed_score) <- c('cell_imputed_score_d0_d10_adaptingFP')
cell_imputed_score$cell_id <- rownames(cell_imputed_score)

comp_df <- merge(fp_d0_d10_treatment, cell_imputed_score, by = 'cell_id')
ggplot(comp_df, aes(x = cell_imputed_score_d0_d10, y = cell_imputed_score_d0_d10_adaptingFP)) +
  geom_point() +
  stat_cor() +
  theme_minimal()
