rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggplot2)
library(ggpubr)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'

remove_unassigned_cells <- TRUE


treatment <- 'CIS'

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

# all_data[['fasttopic_COCL2']] <- all_data_fasttopic_COCL2
# all_data[['peakVI_COCL2']] <- all_data_peakVI_COCL2

all_data[['fasttopic_CIS']] <- all_data_fasttopic_CIS
all_data[['peakVI_CIS']] <- all_data_peakVI_CIS

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
# Calculate the number of week5 cells from fate potential
# =============================================================================
all_data_day10 <- subset(all_data, dataset == paste0('day10_', treatment))
metadat.day10 <- all_data_day10@meta.data
metadat.day10$cell_id <- rownames(metadat.day10)

fp_name <- paste0('fatepotential_', treatment, '_d10_w5')
fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['cell_imputed_score']])
colnames(fp_d10_w5_treatment) <- c('cell_imputed_score')
fp_d10_w5_treatment$cell_id <- rownames(fp_d10_w5_treatment)

fp_d10_w5_treatment$cell_imputed_count <- 10**(fp_d10_w5_treatment$cell_imputed_score)
fp_d10_w5_treatment$cell_imputed_count <- round(fp_d10_w5_treatment$cell_imputed_count, 0)

fp_d10_w5_treatment <- fp_d10_w5_treatment[fp_d10_w5_treatment$cell_imputed_score > 0, ]

fp_d10_w5_treatment <- merge(fp_d10_w5_treatment, metadat.day10[, c('cell_id', 'assigned_lineage')], by = 'cell_id')

imputed_adapting_progeny_size <- fp_d10_w5_treatment %>% 
  group_by(assigned_lineage) %>%
  summarize(n_imputed_adapting_progeny = sum(cell_imputed_count))


acutal_day10_size <- metadat.day10 %>% 
  group_by(assigned_lineage) %>%
  summarize(n_cells = n())


comp_df <- merge(acutal_day10_size, imputed_adapting_progeny_size, by = 'assigned_lineage', all.x = TRUE)

ggplot(comp_df, aes(x = log10(n_cells + 1), y = log10(n_imputed_adapting_progeny + 1))) +
  geom_point() +
  stat_cor(method = 'spearman') +
  labs(x = 'Actual lineage size', y = 'Imputed lineage size') +
  theme_bw()

# =============================================================================
# Calculate day0 lineage size
# =============================================================================
all_data_day0 <- subset(all_data, dataset == 'day0')
metadat.day0 <- all_data_day0@meta.data


input_df <- metadat.day0 %>% 
  group_by(assigned_lineage) %>%
  summarize(current_lineage_size = n())

# =============================================================================
# prepare input for fate potential
# =============================================================================

input_df <- merge(input_df, imputed_adapting_progeny_size, by = 'assigned_lineage', all.x = TRUE)
input_df[is.na(input_df)] <- 0

ggplot(input_df, aes(x = current_lineage_size, y = n_imputed_adapting_progeny)) +
  geom_point() +
  stat_cor() +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'Current lineage size', y = 'Imputed lineage size') +
  theme_bw()

cell_lineage <- all_data_day0$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))

# prepare the input
rna_features <- all_data_day0[[paste0("fasttopic_", treatment)]]@cell.embeddings
atac_features <- all_data_day0[[paste0("peakVI_", treatment)]]@cell.embeddings

cell_features <- cbind(rna_features, atac_features)
cell_features <- scale(cell_features)

# =============================================================================
# Run fate potential
# =============================================================================
set.seed(123)

rownames(input_df) <- input_df$assigned_lineage
input_df <- input_df[, c('current_lineage_size', 'n_imputed_adapting_progeny')]
input_df <- as.matrix(input_df)
colnames(input_df) <- c("now", "future")

lineage_future_count <- input_df[, 'future']

colnames(input_df) <- c("day0", paste0("day10_", treatment))
tab_mat <- input_df

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = paste0("day10_", treatment),
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
  stat_cor(method = 'spearman') +
  # geom_abline(intercept = 0, slope = 1, col = 'red') +
  theme_minimal()

saveRDS(final_fit, '~/Downloads/final_fit_d0_d10_adpating_thres_0_progeny_count_CIS.rds')

