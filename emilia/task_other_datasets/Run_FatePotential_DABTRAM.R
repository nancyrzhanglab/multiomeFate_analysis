rm(list=ls())
library(dplyr)
library(reshape2)
library(Seurat)
library(Signac)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/output/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'Writeup6m_all-data.RData'))

all_data@assays[["geneActivity"]] <- NULL
all_data@assays[["Lineage"]] <- NULL


# ==============================================================================
# Subset data
# ==============================================================================
all_data <- subset(all_data, subset = dataset %in% c('day0', 'day10_DABTRAM'))
all_data[['Is_LB_present']] <- ifelse(is.na(all_data@meta.data$assigned_lineage), 'No', 'Yes')
all_data <- subset(all_data, subset = Is_LB_present == 'Yes')
metadat <- all_data@meta.data
metadat$cell_barcode <- rownames(metadat)

# ==============================================================================
# Prepare for FatePotential
# ==============================================================================
embedding_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings

early_idx <- which(all_data$dataset == 'day0')

lineage_size <- metadat %>% 
  group_by(assigned_lineage, dataset) %>% 
  summarise(num_cells = n())
lineage_size_w <- dcast(lineage_size, assigned_lineage ~ dataset)
lineage_size_w[is.na(lineage_size_w)] <- 0
colnames(lineage_size_w) <- c('assigned_lineage', 'now', 'future')
lineage_size_w <- lineage_size_w[lineage_size_w$now > 0, ]

cell_lineage <- metadat[ ,c('cell_barcode', 'assigned_lineage')]
cell_lineage <- cell_lineage[cell_lineage$cell_barcode %in% names(early_idx), ]

tmp <- cell_lineage$assigned_lineage
names(tmp) <- cell_lineage$cell_barcode
cell_lineage <- tmp
uniq_lineage <- sort(unique(cell_lineage))

lineage_future_count <- lineage_size_w$future
names(lineage_future_count) <- lineage_size_w$assigned_lineage
rownames(lineage_size_w) <- lineage_size_w$assigned_lineage
tab_mat <- lineage_size_w[, c('now', 'future')]

set.seed(10)

fit_res <- multiomeFate:::lineage_cv(
  cell_features = embedding_mat[early_idx,, drop=FALSE],
  cell_lineage = cell_lineage,
  future_timepoint = "future",
  lineage_future_count = lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 10,
  tab_mat = tab_mat,
  num_folds = 10,
  verbose = 2
)

final_fit <- multiomeFate:::lineage_cv_finalize(
  cell_features = embedding_mat[early_idx,,drop=FALSE],
  cell_lineage = cell_lineage,
  fit_res = fit_res,
  lineage_future_count = lineage_future_count
)

lineage_imputed_count_FP <- final_fit$lineage_imputed_count
cell_imputed_score_FP <- final_fit$cell_imputed_score
round(final_fit$coefficient_vec, 5)
cell_imputed_count_FP <- 10**(cell_imputed_score)
hist(cell_imputed_score_FP, breaks = 100, main = 'Histogram of 10^(cell_imputed_score)')

# ==============================================================================
# Compare with Kevin's output
# ==============================================================================
load(paste0(in_dir, 'Growth_potential/Writeup6r_DABTRAM_day0_lineage-imputation_postprocess.RData'))

cell_imputed_score_df <- as.data.frame(cell_imputed_score)
cell_imputed_score_FP_df <- as.data.frame(cell_imputed_score_FP)

cell_imputed_score_df$cell_id <- rownames(cell_imputed_score_df)
cell_imputed_score_FP_df$cell_id <- rownames(cell_imputed_score_FP_df)

comp <- merge(cell_imputed_score_df, cell_imputed_score_FP_df, by = 'cell_id')

ggplot(comp, aes(x = cell_imputed_score, y = cell_imputed_score_FP)) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  xlab('Writeup6r_DABTRAM_day0_lineage-imputation_postprocess') +
  ylab('FP') +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5)) +
  theme_bw()

