rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)

data_dir <- '/Users/emiliac/Library/CloudStorage/Dropbox/Thesis/Lineage_trace/data_other/Weinreb_Larry_Hematopoiesis/'
output_dir <- '/Users/emiliac/Library/CloudStorage/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/Weinreb_Larry_Hematopoiesis/'
date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

SAMPLE_NAME <- 'Weinreb_Larry_Hematopoiesis'
DAY <- 'd6'
# =============================================================================
# Load data
# =============================================================================
load(paste0(data_dir, "/Writeup13_combined_fatepotential.RData"))

# =============================================================================
# Wrangle data
# =============================================================================
metadat <- seurat_obj@meta.data

# grep a pattern from the library name
metadat$day <- gsub('.*d([0-9]+).*', '\\1', metadat$Library)
metadat$day <- paste0('d', metadat$day)

seurat_obj <- AddMetaData(seurat_obj, metadat)

# subset to day
seurat_obj.use <- subset(seurat_obj, subset = day == DAY)
metadat.use <- seurat_obj.use@meta.data
metadat.use$cell_barcode <- rownames(metadat.use)

# get PCA
data_use <- seurat_obj.use@reductions[["pca"]]@cell.embeddings
data_use <- as.data.frame(data_use)

NUM_DIM = dim(data_use)[2]
MIN_CELL = 4

# merge with metadata
data_use$cell_barcode <- rownames(data_use)

data_use <- merge(data_use, metadat.use[, c('cell_barcode', 'assigned_lineage')], by = "cell_barcode")

# calculate lineage size
lineage_barcodes = data_use %>% 
  group_by(assigned_lineage) %>% 
  summarise(n_cells = n()) %>% 
  arrange(desc(n_cells))

# =============================================================================
# Calculate pairwise distance within a lineage barcode 
# =============================================================================

d_in_df = data.frame(matrix(nrow = 0, ncol = 3))
colnames(d_in_df) = c('assigned_lineage', 'avg_euc_dist', 'n_cells')

embedding_col_names = paste0('PC_', 1:NUM_DIM)

for(l in unique(data_use$assigned_lineage)) {
  df = data_use %>% filter(assigned_lineage == l)
  if(nrow(df) >= MIN_CELL) {
    dists = dist(df[, embedding_col_names], method = "euclidean")
    avg_euc_dist = mean(dists)
  }else {
    avg_euc_dist = 0
  }
  d_in_df = rbind(d_in_df, 
                  data.frame(assigned_lineage = l, 
                             avg_euc_dist = avg_euc_dist, 
                             n_cells = nrow(df)))
  
}

d_in_df = d_in_df[d_in_df$n_cells >= MIN_CELL, ]

# =============================================================================
# Calculate pairwise distance shuffled barcodes
# =============================================================================
data_use$shuffled_assigned_lineage = sample(data_use$assigned_lineage, size = nrow(data_use), replace = F)

d_in_shuffled_df = data.frame(matrix(nrow = 0, ncol = 3))
colnames(d_in_shuffled_df) = c('assigned_lineage', 'avg_euc_dist', 'n_cells')

for(l in unique(data_use$assigned_lineage)) {
  df = data_use %>% filter(shuffled_assigned_lineage == l)
  if(nrow(df) >= MIN_CELL) {
    dists = dist(df[, embedding_col_names], method = "euclidean")
    avg_euc_dist = mean(dists)
  }else {
    avg_euc_dist = 0
  }
  d_in_shuffled_df = rbind(d_in_shuffled_df, 
                           data.frame(assigned_lineage = l, 
                                      avg_euc_dist = avg_euc_dist, 
                                      n_cells = nrow(df)))
}

d_in_shuffled_df = d_in_shuffled_df[d_in_shuffled_df$n_cells >= MIN_CELL, ]

# =============================================================================
# Process Results
# =============================================================================
print(dim(d_in_df))
print(dim(d_in_shuffled_df))

d_in_df$normalized_avg_eud_dist_by_shuffle = d_in_df$avg_euc_dist / median(d_in_shuffled_df$avg_euc_dist)
d_in_shuffled_df$normalized_avg_eud_dist_by_shuffle = d_in_shuffled_df$avg_euc_dist / median(d_in_shuffled_df$avg_euc_dist)

write.csv(d_in_df, paste0(output_dir, '/lineage_variability_', SAMPLE_NAME, '_', DAY, '.csv'), row.names = F)
write.csv(d_in_shuffled_df, paste0(output_dir, '/lineage_variability_shuffled_', SAMPLE_NAME,'_', DAY, '.csv'), row.names = F)

# =============================================================================
# Plotting
# =============================================================================

# Plotting lineage size distribution
p1 = ggplot(lineage_barcodes, aes(x = n_cells)) + 
  geom_histogram(binwidth = 5) + 
  theme_bw() +
  labs(title = paste0('Lineage size distribution (', SAMPLE_NAME, ', ', DAY, ')'), y = 'Number of cells')
ggsave(paste0(output_dir, '/lineage_size_histogram_', SAMPLE_NAME, '_', DAY, '.png'), plot = p1, width = 5, height = 4 )

# Plotting lineage variability in histogram
d_in_df$category = 'Data'
d_in_shuffled_df$category = 'Shuffled'

d_in_plot = rbind(d_in_df, d_in_shuffled_df)
p2 = ggplot(d_in_plot, aes(x = normalized_avg_eud_dist_by_shuffle)) + 
  geom_histogram(aes(y = ..density.., fill = category), col = 'black', alpha = 0.5, position = 'identity') +
  geom_density(aes(col = category), lwd = 1) +
  # xlim(-0.1, 3) +
  # coord_cartesian(xlim = c(0.2, 1.2)) +
  theme_bw() +
  labs(title = paste0('Lineage variability (', SAMPLE_NAME, ', ', DAY, ')'), 
       x = 'Normalized average Euclidean distance')
ggsave(paste0(output_dir, '/normalized_avg_eud_dist_', SAMPLE_NAME, '_', DAY, '.png'), plot = p2, width = 6, height = 4 )

# Plotting lineage variability's correlation with lineage size
p3 = ggplot(d_in_df, aes(x = n_cells, y = normalized_avg_eud_dist_by_shuffle)) + 
  geom_point() + 
  theme_bw() +
  labs(title = paste0('Lineage variability vs lineage size (', SAMPLE_NAME, ', ', DAY, ')'), 
       x = 'Number of cells', y = 'Normalized average Euclidean distance')
ggsave(paste0(output_dir, '/lineage_variability_vs_size_', SAMPLE_NAME, '_', DAY, '.png'), plot = p3, width = 5, height = 5 )

