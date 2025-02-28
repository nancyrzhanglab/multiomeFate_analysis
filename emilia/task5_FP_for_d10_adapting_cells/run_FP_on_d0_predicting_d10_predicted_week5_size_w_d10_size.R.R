rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)

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
# Get predicted week5 lineage size
# =============================================================================

all_data_day10 <- subset(all_data, dataset == paste0('day10_', treatment))
metadat.day10 <- all_data_day10@meta.data
metadat.day10$cell_id <- rownames(metadat.day10)

fp_name <- paste0('fatepotential_', treatment, '_d10_w5')

fp_d10_w5_treatment <- as.data.frame(all_data_day10@misc[[fp_name]][['lineage_imputed_count']])
colnames(fp_d10_w5_treatment) <- c('lineage_imputed_count')
fp_d10_w5_treatment$assigned_lineage <- rownames(fp_d10_w5_treatment)

quantile(fp_d10_w5_treatment$lineage_imputed_count)
hist(fp_d10_w5_treatment$lineage_imputed_count, breaks = 50)

# =============================================================================
# Get day10 lineage size
# =============================================================================
lineage_size <- metadat.day10 %>% 
  group_by(assigned_lineage) %>%
  summarize(n_cells.day10 = n())


# =============================================================================
# Prepare input for fate potential
# =============================================================================

# prepare the response variables
all_data_day0 <- subset(all_data, dataset == 'day0')
metadat.day0 <- all_data_day0@meta.data
metadat.day0$cell_id <- rownames(metadat.day0)

input_df <- metadat.day0 %>% 
  group_by(assigned_lineage) %>%
  summarize(lineage_current_count = n())

lineage_size.day10 <- merge(metadat.day0[, c('cell_id', 'assigned_lineage')], lineage_size, by = 'assigned_lineage', all.x = TRUE)
lineage_size.day10[is.na(lineage_size.day10)] <- 0
rownames(lineage_size.day10) <- lineage_size.day10$cell_id


input_df <- merge(input_df, fp_d10_w5_treatment, by = 'assigned_lineage', all.x = TRUE)
input_df[is.na(input_df)] <- 0


ggplot(input_df, aes(x = lineage_current_count, y = lineage_imputed_count)) +
  geom_jitter(width = 0.1) +
  stat_cor() +
  theme_minimal()


cell_lineage <- all_data_day0$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))

# prepare the input
rna_features <- all_data_day0[[paste0("fasttopic_", treatment)]]@cell.embeddings
atac_features <- all_data_day0[[paste0("peakVI_", treatment)]]@cell.embeddings
lineage_size.day10 <- lineage_size.day10[rownames(atac_features), ]
lineage_size.day10 <- lineage_size.day10$n_cells.day10
names(lineage_size.day10) <- rownames(atac_features)

cell_features <- cbind(rna_features, atac_features, lineage_size.day10)
cell_features <- scale(cell_features)
# =============================================================================
# Run fate potential
# =============================================================================
set.seed(123)

rownames(input_df) <- input_df$assigned_lineage
input_df <- input_df[, c('lineage_current_count', 'lineage_imputed_count')]
input_df <- as.matrix(input_df)
colnames(input_df) <- c("now", "future")

lineage_future_count <- input_df[, 'future']

colnames(input_df) <- c("day0", paste0("day10_imputed_w5_", treatment))
tab_mat <- input_df

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = paste0("day10_imputed_w5_", treatment),
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
ggplot(comp_df, aes(x = log10(lineage_future_count+1), y = log10(lineage_imputed_count + 1))) +
  geom_jitter(width = 0.05) +
  stat_cor() +
  # geom_abline(intercept = 0, slope = 1, col = 'red') +
  theme_minimal()

View(cell_imputed_score)

saveRDS(final_fit, '~/Downloads/final_fit_d0_d10_predicted_w5_w_d10_size.CIS_1.rds')

cell_features2 <- cell_features[, -c(ncol(cell_features))]
cell_features2 <- cbind(1, cell_features2)
colnames(cell_features2)[1] <- "Intercept"
final_fit.coefficient_vec2 <- as.matrix(final_fit$coefficient_vec[-c(length(final_fit$coefficient_vec))])

cell_imputed_score2 <- as.numeric(cell_features2 %*% final_fit.coefficient_vec2)
hist(cell_imputed_score2)
names(cell_imputed_score2) <- rownames(cell_features)

saveRDS(cell_imputed_score2, '~/Downloads/cell_imputed_score2_d0_d10_predicted_w5_w_d10_size.CIS_1.rds')
