rm(list=ls())
library(Seurat)
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

in_dir <- '~/Dropbox/Thesis/Lineage_trace/data/Watermelon/'
out_dir <- '~/Dropbox/Thesis/Lineage_trace/output/Watermelon/'

# ==============================================================================
# Read data
# ==============================================================================
load(paste0(in_dir, 'PC9_time_course_fasttopics.RData'))
# Seurat::Idents(seurat_object) <- "time_point"

# ==============================================================================
# Subset data
# ==============================================================================
seurat_subset <- subset(seurat_object, subset = time_point %in% c(7, 14))
seurat_subset[['Is_LB_present']] <- ifelse(is.na(seurat_subset@meta.data$lineage_barcode), 'No', 'Yes')
seurat_subset <- subset(seurat_subset, subset = Is_LB_present == 'Yes')
metadat <- seurat_subset@meta.data
metadat$cell_barcode <- rownames(metadat)
metadat$cell_ID <- paste0('Day', metadat$time_point, '-', metadat$cell_barcode)

# ==============================================================================
# Prepare for FatePotential
# ==============================================================================
embedding_mat <- seurat_subset[["fasttopic"]]@cell.embeddings

early_idx <- which(seurat_subset$time_point == 7)

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


plot(x = log10(lineage_size_w$future + 1),
     y = log10(lineage_imputed_count[lineage_size_w$lineage_barcode] + 1), 
     asp = TRUE, pch = 16,
     xlab = "Observed lineage size (Day 14)",
     ylab = "Fitted lineage size (from Day 7 to Day 14)",
     main = paste0("Correlation: ", round(stats::cor(
       lineage_size_w$future,
       lineage_imputed_count[lineage_size_w$lineage_barcode]
     ), 2))
)


metadat <- seurat_object@meta.data
metadat$cell_barcode <- rownames(metadat)
final_fit_df <- as.data.frame(final_fit[["cell_imputed_score"]])
colnames(final_fit_df) <- 'FatePotential'
final_fit_df$cell_barcode <- rownames(final_fit_df)
final_fit_df$time_point <- 7

metadat <- merge(metadat, final_fit_df, by = c('cell_barcode', 'time_point'), all = T)
metadat$FatePotential[is.na(metadat$FatePotential)] <- 0
seurat_object <- AddMetaData(seurat_object, metadat)

FeaturePlot(seurat_object, 'FatePotential')

FeaturePlot(seurat_object, 'time_point')

umap_coord <- as.data.frame(seurat_object@reductions[["umap"]]@cell.embeddings)
umap_coord$cell_id <- rownames(umap_coord)

metadat <- merge(metadat, umap_coord, by.x = 'cell_barcode', by.y = 'cell_id')
metadat_7 <- metadat[metadat$time_point == 7, ]
ggplot() +
  geom_point(data = metadat, aes(x = umap_1, y = umap_2), color = '#F6DCAC', alpha=0.5, size = 0.5) +
  geom_point(data = metadat_7, aes(x = umap_1, y = umap_2, color = FatePotential), size = 0.5) +
  scale_color_gradient(low = 'blue', high = 'red') +
  theme_classic()
