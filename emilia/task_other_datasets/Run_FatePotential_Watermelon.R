rm(list=ls())
library(dplyr)
library(reshape2)
library(Seurat)
library(Signac)
library(multiomeFate)
library(ggplot2)
library(ggpubr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Watermelon/'

now_timepoint <- 3
future_timepoint <- 14

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
# Seurat::Idents(seurat_object) <- "time_point"

# ==============================================================================
# Subset data
# ==============================================================================
seurat_subset <- subset(seurat_object, subset = time_point %in% c(now_timepoint, future_timepoint))
seurat_subset[['Is_LB_present']] <- ifelse(is.na(seurat_subset@meta.data$lineage_barcode), 'No', 'Yes')
seurat_subset <- subset(seurat_subset, subset = Is_LB_present == 'Yes')
metadat <- seurat_subset@meta.data
metadat$cell_barcode <- rownames(metadat)
metadat$cell_ID <- paste0('Day', metadat$time_point, '-', metadat$cell_barcode)

# ==============================================================================
# Prepare for FatePotential
# ==============================================================================
embedding_mat <- seurat_subset[["fasttopic"]]@cell.embeddings
embedding_mat <- scale(embedding_mat)

early_idx <- which(seurat_subset$time_point == now_timepoint)

lineage_size <- metadat %>% 
  group_by(lineage_barcode, time_point) %>% 
  summarise(num_cells = n())
lineage_size_w <- dcast(lineage_size, lineage_barcode ~ time_point)
lineage_size_w[is.na(lineage_size_w)] <- 0
colnames(lineage_size_w) <- c('lineage_barcode', 'now', 'future')
lineage_size_w <- lineage_size_w[lineage_size_w$now > 0, ]

cell_lineage <- metadat[ ,c('cell_barcode', 'lineage_barcode')]
cell_lineage <- cell_lineage[cell_lineage$cell_barcode %in% names(early_idx), ]

tmp <- cell_lineage$lineage_barcode
names(tmp) <- cell_lineage$cell_barcode
cell_lineage <- tmp
uniq_lineage <- sort(unique(cell_lineage))

lineage_future_count <- lineage_size_w$future
names(lineage_future_count) <- lineage_size_w$lineage_barcode
rownames(lineage_size_w) <- lineage_size_w$lineage_barcode
tab_mat <- lineage_size_w[, c('now', 'future')]

# ggplot(tab_mat, aes(x = now, y = future)) +
#   geom_point()

set.seed(10)
fit_res <- multiomeFate:::lineage_cv(
  cell_features = embedding_mat[early_idx,, drop=FALSE],
  cell_lineage = cell_lineage,
  future_timepoint = "future",
  lineage_future_count = lineage_future_count,
  lambda_initial = 5,
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

lineage_imputed_count <- final_fit$lineage_imputed_count
cell_imputed_score <- final_fit$cell_imputed_score
round(final_fit$coefficient_vec, 5)
cell_imputed_count <- 10**(cell_imputed_score)
hist(cell_imputed_count, breaks = 100, main = 'Histogram of 10^(cell_imputed_score)')


plot(x = log10(lineage_size_w$future + 1),
     y = log10(lineage_imputed_count[lineage_size_w$lineage_barcode] + 1), 
     asp = TRUE, pch = 16,
     xlab = "Observed lineage size (Day 14)",
     ylab = paste0("Fitted lineage size (from Day", now_timepoint, " to Day ", future_timepoint, " )"),
     main = paste0("Correlation: ", round(stats::cor(
       lineage_size_w$future,
       lineage_imputed_count[lineage_size_w$lineage_barcode]
     ), 2))
)

lineage_imputed_count_df <- as.data.frame(lineage_imputed_count)
lineage_imputed_count_df$lineage_barcode <- rownames(lineage_imputed_count_df)
lineage_imputed_count_df <- merge(lineage_imputed_count_df, lineage_size_w, by = 'lineage_barcode')
ggplot(lineage_imputed_count_df, 
       aes( x = log10(future + 1), 
            y = log10(lineage_imputed_count + 1))) +
  geom_jitter(alpha = 0.5, width = 0.02) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)+
  ylab(paste0('Log10(imputed lineage size (from Day ', now_timepoint, 'to Day ', future_timepoint, ' )')) +
  xlab(paste0('Log10(observed Day ', future_timepoint, ' count )')) +
  theme_bw() +
  theme(panel.grid = element_blank())


metadat <- seurat_object@meta.data
# metadat <- subset(metadat, select = -c(FatePotential))
metadat$cell_barcode <- rownames(metadat)
final_fit_df <- as.data.frame(final_fit[["cell_imputed_score"]])
final_fit_df <- as.data.frame(cell_imputed_count)
colnames(final_fit_df) <- 'FatePotential'
final_fit_df$cell_barcode <- rownames(final_fit_df)
final_fit_df$time_point <- now_timepoint

metadat <- merge(metadat, final_fit_df, by = c('cell_barcode', 'time_point'), all = T)
metadat$FatePotential[is.na(metadat$FatePotential)] <- 0
# metadat$FatePotential[is.na(metadat$cell_imputed_count)] <- -1
seurat_object <- AddMetaData(seurat_object, metadat)

# FeaturePlot(seurat_object, 'FatePotential')

umap_coord <- as.data.frame(seurat_object@reductions[["umap"]]@cell.embeddings)
umap_coord$cell_id <- rownames(umap_coord)

metadat <- merge(metadat, umap_coord, by.x = 'cell_barcode', by.y = 'cell_id')
metadat_now <- metadat[metadat$time_point == now_timepoint, ]
ggplot() +
  geom_point(data = metadat, aes(x = umap_1, y = umap_2), color = '#F6DCAC', alpha=0.5, size = 0.5) +
  geom_point(data = metadat_now, aes(x = umap_1, y = umap_2, color = FatePotential), size = 0.5) +
  scale_color_gradient(low = 'blue', high = 'red') +
  # xlim(3, 10) +
  # ylim(-5, 5) +
  theme_classic()


save(final_fit, date_of_run, 
     file = paste0(out_dir, 'PC9_time_course_d', now_timepoint, '_d', future_timepoint, '.RData'))

