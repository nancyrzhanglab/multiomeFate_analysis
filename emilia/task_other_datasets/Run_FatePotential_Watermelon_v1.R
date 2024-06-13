rm(list=ls())
library(Seurat)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/outputs/task5_cospar/Watermelon/'

d14_fate <- 'high' # high = non-cycling, median = mid, low = cycling
current <- 0
future <- 14

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))

# ==============================================================================
# Subset data
# ==============================================================================

# Get D14 index
d14_idx <- intersect(which(seurat_object$time_point == future),
                     which(seurat_object$sample_type != paste0('14_', d14_fate)))
d14_idx <- intersect(d14_idx, which(!is.na(seurat_object$lineage_barcode)))

# Get D7 index
d7_idx <- intersect(which(seurat_object$time_point == current), which(!is.na(seurat_object$lineage_barcode)))

# Combine and subset main dataset
keep_idx <- c(d14_idx, d7_idx)
keep_vec <- rep(FALSE, ncol(seurat_object))
keep_vec[keep_idx] <- TRUE
seurat_object$keep <- keep_vec

seurat_subset <- subset(seurat_object, keep == TRUE)
metadat <- seurat_subset@meta.data
metadat$cell_barcode <- rownames(metadat)
metadat$cell_ID <- paste0('Day', metadat$time_point, '-', metadat$cell_barcode)

print(unique(metadat$time_point))
print(unique(metadat$sample_name))
print(dim(metadat))

metadat.all <- seurat_object@meta.data
metadat.all <- metadat.all[metadat.all$time_point %in% c(current, future), ]
metadat.all <- metadat.all[metadat.all$sample_type != '14_high', ]
metadat.all <- metadat.all %>% drop_na()

print(unique(metadat.all$time_point))
print(unique(metadat.all$sample_name))
print(dim(metadat.all))


# ==============================================================================
# Prepare for FatePotential
# ==============================================================================
embedding_mat <- seurat_subset[["fasttopic"]]@cell.embeddings

# par(mfrow = c(4, 2))
# hist(embedding_mat[, 1], breaks = 100, main = 'FastTopic_1_COCL2')
# hist(embedding_mat_1[, 1], breaks = 100, main = 'FastTopic_1_PC9')
# hist(embedding_mat[, 2], breaks = 100, main = 'FastTopic_2_COCL2')
# hist(embedding_mat_1[, 2], breaks = 100, main = 'FastTopic_2_PC9')
# hist(embedding_mat[, 3], breaks = 100, main = 'FastTopic_3_COCL2')
# hist(embedding_mat_1[, 3], breaks = 100, main = 'FastTopic_3_PC9')
# hist(embedding_mat[, 4], breaks = 100, main = 'FastTopic_4_COCL2')
# hist(embedding_mat_1[, 4], breaks = 100, main = 'FastTopic_4_PC9')

early_idx <- which(seurat_subset$time_point == current)

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

lineage_imputed_count <- final_fit$lineage_imputed_count
cell_imputed_score <- final_fit$cell_imputed_score
round(final_fit$coefficient_vec, 5)
cell_imputed_count <- 10**(cell_imputed_score)
hist(cell_imputed_count, breaks = 100, main = 'Histogram of 10^(cell_imputed_score)')

plot(x = log10(lineage_size_w$future + 1),
     y = log10(lineage_imputed_count[lineage_size_w$lineage_barcode] + 1), 
     asp = TRUE, pch = 16,
     xlab = "Observed lineage size (Day 14 mid / low)",
     ylab = "Fitted lineage size (from Day 0 to Day 14 mid / low)",
     main = paste0("Correlation: ", round(stats::cor(
       lineage_size_w$future,
       lineage_imputed_count[lineage_size_w$lineage_barcode]
     ), 2))
)

metadat <- seurat_object@meta.data
metadat$cell_barcode <- rownames(metadat)
# final_fit_df <- as.data.frame(final_fit[["cell_imputed_score"]])
final_fit_df <- as.data.frame(cell_imputed_count)
# colnames(final_fit_df) <- 'FatePotential'
final_fit_df$cell_barcode <- rownames(final_fit_df)
final_fit_df$time_point <- current

metadat <- merge(metadat, final_fit_df, by = c('cell_barcode', 'time_point'), all = T)
metadat$cell_imputed_count[is.na(metadat$cell_imputed_count)] <- 0
seurat_object <- AddMetaData(seurat_object, metadat)

umap_coord <- as.data.frame(seurat_object@reductions[["umap"]]@cell.embeddings)
umap_coord$cell_id <- rownames(umap_coord)

metadat <- merge(metadat, umap_coord, by.x = 'cell_barcode', by.y = 'cell_id')
metadat_now <- metadat[metadat$time_point == 0, ]
ggplot() +
  geom_point(data = metadat, aes(x = umap_1, y = umap_2), color = '#F6DCAC', alpha=0.5, size = 0.5) +
  geom_point(data = metadat_now, aes(x = umap_1, y = umap_2, color = cell_imputed_count), size = 0.5) +
  scale_color_gradient(low = 'blue', high = 'red') +
  # xlim(3, 10) +
  # ylim(-5, 5) +
  theme_classic()

d14_fate <- 'mid_low'
save(final_fit, date_of_run, session_info, file = paste0(out_dir, 'PC9_time_course_d', current, '_to_d14_', d14_fate, '.RData'))

