rm(list=ls())
library(Seurat)
library(multiomeFate)
library(tidyverse)
library(ggpubr)

# =============================================================================
# reading data
# =============================================================================

# For in vivo data
data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Minn_lab/Final/'
out_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task7_run_FP_on_R499_TBK1i_Lineage_Trace/'
base_name <- 'TBK1i_Multiome_StringentMT_InVivo_Run3_ATAC_Cancer_Harmony_eigs18'

metadat <- readRDS(paste0(data_dir, base_name, '/', base_name, '.scMultiome_combined_object_METADAT.rds'))
snapATAC.dimrec <- read.csv(paste0(data_dir, base_name, '/', 'Dimensionality_Reduction_Harmony.', base_name, '.csv'), row.names = 1)

# check if all barcodes the same
all(rownames(snapATAC.dimrec) == rownames(metadat))

# =============================================================================
# Prepare the response variables
# =============================================================================
metadat <- metadat[!is.na(metadat$lineage_barcode_assignment), ]
metadat <- metadat[metadat$lineage_barcode_assignment != 'Too many barcodes', ]
metadat.JAKi_TBK1i.dICB <- metadat[metadat$Condition == 'JAKi TBK1i dICB D19 InVivo', ]

lineage_size <- metadat.JAKi_TBK1i.dICB %>% 
  group_by(lineage_barcode_assignment) %>%
  summarize(n_cells = n())

quantile(lineage_size$n_cells)
hist(lineage_size$n_cells, breaks = 50)

# =============================================================================
# Prepare input for fate potential
# =============================================================================

metadat.JAKi_TBK1i.d14 <- metadat[metadat$Condition == 'JAKi TBK1i D14 InVivo', ]
snapATAC.dimrec.JAKi_TBK1i.d14 <- snapATAC.dimrec[rownames(metadat.JAKi_TBK1i.d14), ]

input_df <- metadat.JAKi_TBK1i.d14 %>% 
  group_by(lineage_barcode_assignment) %>%
  summarize(lineage_current_count = n())


input_df <- merge(input_df, lineage_size, by = 'lineage_barcode_assignment', all.x = TRUE)
input_df[is.na(input_df)] <- 0


ggplot(input_df, aes(x = log10(lineage_current_count+1), y = log10(n_cells+1))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = 'red') +
  stat_cor() +
  theme_minimal()

cell_lineage <- metadat.JAKi_TBK1i.d14$lineage_barcode_assignment
uniq_lineage <- sort(unique(cell_lineage))

metadat.JAKi_TBK1i.d14$row.names <- rownames(metadat.JAKi_TBK1i.d14)
metadat.JAKi_TBK1i.d14.cur.size <- metadat.JAKi_TBK1i.d14[, c('row.names', 'lineage_barcode_assignment')]
metadat.JAKi_TBK1i.d14.cur.size <- merge(metadat.JAKi_TBK1i.d14.cur.size, input_df, by = 'lineage_barcode_assignment')
rownames(metadat.JAKi_TBK1i.d14.cur.size) <- metadat.JAKi_TBK1i.d14.cur.size$row.names
metadat.JAKi_TBK1i.d14.cur.size <- metadat.JAKi_TBK1i.d14.cur.size[rownames(metadat.JAKi_TBK1i.d14), ]

# check if all
all(rownames(snapATAC.dimrec.JAKi_TBK1i.d14) == rownames(metadat.JAKi_TBK1i.d14.cur.size))

# prepare the input
snapATAC.dimrec.JAKi_TBK1i.d14 <- snapATAC.dimrec[rownames(metadat.JAKi_TBK1i.d14), ]
metadat.JAKi_TBK1i.d14.cur.size <- metadat.JAKi_TBK1i.d14.cur.size %>% pull(lineage_current_count)
# cell_features <- as.matrix(cbind(snapATAC.dimrec.JAKi_TBK1i.d14, metadat.JAKi_TBK1i.d14.cur.size))
cell_features <- as.matrix(snapATAC.dimrec.JAKi_TBK1i.d14)
# cell_features <- scale(cell_features)

# =============================================================================
# Run fate potential
# =============================================================================
set.seed(123)

rownames(input_df) <- input_df$lineage_barcode_assignment
input_df <- input_df[, c('lineage_current_count', 'n_cells')]
input_df <- as.matrix(input_df)
colnames(input_df) <- c("now", "future")

lineage_future_count <- input_df[, 'future']

colnames(input_df) <- c("day0", 'dICB_JAKi_TBK1i')
tab_mat <- input_df

fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = 'dICB_JAKi_TBK1i',
  lineage_future_count = lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 20,
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
  geom_jitter(width = 0.1) +
  stat_cor() +
  # geom_abline(intercept = 0, slope = 1, col = 'red') +
  theme_minimal()

View(cell_imputed_score)

saveRDS(final_fit, '~/Downloads/final_fit_d14_d19_JAKi_TBK1i.rds')

cell_features2 <- cell_features[, -c(ncol(cell_features))]
cell_features2 <- cbind(1, cell_features2)
colnames(cell_features2)[1] <- "Intercept"
final_fit.coefficient_vec2 <- as.matrix(final_fit$coefficient_vec[-c(length(final_fit$coefficient_vec))])

cell_imputed_score2 <- as.numeric(cell_features2 %*% final_fit.coefficient_vec2)
hist(cell_imputed_score2)
names(cell_imputed_score2) <- rownames(cell_features)

saveRDS(cell_imputed_score2, '~/Downloads/final_fit_d14_d19_JAKi_TBK1i_w_d14_size_scaled.fatepotential.rds')

# =============================================================================
# Compare
# =============================================================================
input_df <- metadat.JAKi_TBK1i.d14 %>% 
  group_by(lineage_barcode_assignment) %>%
  summarize(lineage_current_count = n())
input_df <- merge(input_df, lineage_size, by = 'lineage_barcode_assignment', all.x = TRUE)
input_df[is.na(input_df)] <- 0
colnames(input_df) <- c('lineage_barcode_assignment', 'lineage_current_count', 'lineage_future_count')

comp_df2 <- merge(comp_df, input_df, by.x = 'assigned_lineage', by.y = 'lineage_barcode_assignment')

ggplot(comp_df2, aes(x = log10(lineage_future_count.x + 1), y = log10(lineage_imputed_count + 1)))+
  geom_jitter(width = 0.1) +
  stat_cor() +
  theme_minimal()

write.csv(comp_df2, paste0(out_dir, 'JAKi_TBK1i_d14_d19_lineage_size_comparison.csv'), row.names = F)

