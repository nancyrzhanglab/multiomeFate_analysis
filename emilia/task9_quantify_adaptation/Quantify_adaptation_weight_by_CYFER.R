rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
score_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'


remove_unassigned_cells <- TRUE

treatment <- 'CIS'
time1 <- 'day10'
time2 <- 'week5'

# =============================================================================
# Read data
# =============================================================================
load(paste0(data_dir, 'Writeup10a_data_empty.RData'))
load(paste0(data_dir, 'Writeup10a_data_rna_dimred.RData'))
load(paste0(data_dir, 'Writeup10a_data_fatepotential.RData'))
load(paste0(data_dir, 'Writeup10a_data_fasttopic_', treatment, '.RData'))

all_data[["pca"]] <- all_data_pca
all_data[["umap"]] <- all_data_umap
all_data@misc <- all_data_fatepotential
# all_data[['fasttopic_DABTRAM']] <- all_data_fasttopic_DABTRAM
# all_data[['fasttopic_COCL2']] <- all_data_fasttopic_COCL2

# remove cells with no lineage
if(remove_unassigned_cells) {
  print("Removing cells with no assigned lineage")
  all_data$keep <- !is.na(all_data$assigned_lineage)
  if(any(!all_data$keep)){
    print(paste0("There are ", length(which(!all_data$keep)), " cells being removed"))
    all_data <- subset(all_data, keep == TRUE)
  }
}

metadat <- all_data@meta.data
umap <- all_data@reductions[["umap"]]@cell.embeddings
pca <- all_data@reductions[["pca"]]@cell.embeddings
# fasttopic <- all_data_fasttopic_DABTRAM@cell.embeddings

# =============================================================================
# Subset data
# =============================================================================

# subset metadat by time point
if(time1 != 'day0') {
  metadat.time1 <- metadat[metadat$dataset == paste0(time1, '_', treatment), ]
}else {
  metadat.time1 <- metadat[metadat$dataset == 'day0', ]
}

metadat.time2 <- metadat[metadat$dataset == paste0(time2, '_', treatment), ]  

lin.common <- intersect(metadat.time1$assigned_lineage, metadat.time2$assigned_lineage)

# subset metadat by lineage
metadat.time1 <- metadat.time1[metadat.time1$assigned_lineage %in% lin.common, ]
metadat.time2 <- metadat.time2[metadat.time2$assigned_lineage %in% lin.common, ]

pca.time1 <- pca[rownames(pca) %in% rownames(metadat.time1), ] # fasttopic[rownames(fasttopic) %in% rownames(metadat.time1), ]
pca.time2 <- pca[rownames(pca) %in% rownames(metadat.time2), ] # fasttopic[rownames(fasttopic) %in% rownames(metadat.time2), ] 

# =============================================================================
# Helper to calculate pairwise euclidean distance
# =============================================================================
euclidean_dist <- function(dimrec_t1, dimrec_t2) {
  # for each row in dimrec_t1, calculate the euclidean distance to each row in dimrec_t2
  # and return a matrix of distances
  
  # Ensure A and B are matrices
  A <- as.matrix(dimrec_t1)
  B <- as.matrix(dimrec_t2)
  
  dist_sq <- apply(A, MARGIN=1, function(a) {
    # for each row in B, subtract corresponding entry in a
    # and square the result
    aa <- matrix(rep(a, nrow(B)), nrow=nrow(B), byrow=TRUE)
    b_a <- B - aa
    b_a_sq <- b_a^2
    # sum the squared differences
    b_a_sq_sum <- rowSums(b_a_sq)
    })
  
  
  # Return the Euclidean distances
  sqrt(dist_sq)
}


weighted_mean_by_fatepotential <- function(dist_mat.percell.mean, fatepotential_t1_fot_t2.lin.df) {
  
  # make name consistant
  dist_mat.percell.mean <- as.data.frame(dist_mat.percell.mean)
  dist_mat.percell.mean$cell_id <- rownames(dist_mat.percell.mean)
  
  fatepotential_t1_fot_t2.lin.df$cell_id <- rownames(fatepotential_t1_fot_t2.lin.df)
  
  # make dataframe with dist_mat.percell.mean, fatepotential
  df <- merge(dist_mat.percell.mean, fatepotential_t1_fot_t2.lin.df, by = 'cell_id')
  colnames(df) <- c('cell_id', 'dist_mat.percell.mean', 'fatepotential')
  df$fatepotential.shift <- as.numeric(df$fatepotential) - min(df$fatepotential)
  df$progeny_count <- 10**df$fatepotential
  
  # calcualte fatepotential weight
  df$fp_weight <- df$progeny_count / sum(df$progeny_count)
  # df$fp_weight <- df$fatepotential.shift / sum(df$fatepotential.shift)
  
  # weight distance by fatepotential
  df$fp_weighted_dist <- df$fp_weight * df$dist_mat.percell.mean
  
  return(sum(df$fp_weighted_dist))
}

match_name <- function(time){
  time <- gsub('day', 'd', time)
  time <- gsub('week', 'w', time)
}

# =============================================================================
# Distance between time points
# =============================================================================

# within lineages
dist_df <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(dist_df) <- c('lineage', 'dist', 'n_cells.t1', 'n_cells.t2')

for (lin in lin.common) {
  print(lin)
  cell.time1 <- rownames(metadat.time1[metadat.time1$assigned_lineage == lin, ])
  cell.time2 <- rownames(metadat.time2[metadat.time2$assigned_lineage == lin, ])
  
  pca.time1.lin <- pca.time1[cell.time1, ]
  pca.time2.lin <- pca.time2[cell.time2, ]
  
  if(length(cell.time1) < 2 | length(cell.time2) < 2) {
    print(paste0("Not enough cells in lineage ", lin, " at time point ", time1, " or ", time2))
    next
  }
  
  # Calculate the pairwise distances
  dist_mat <- euclidean_dist(pca.time1.lin, pca.time2.lin)
  
  # Calculate row mean
  dist_mat.percell.mean <- colMeans(dist_mat, na.rm = TRUE)
  
  # Get the mean distance for this lineage
  fatepotential_t1_fot_t2.lin.df <- metadat.time1[metadat.time1$assigned_lineage == lin, ] %>% 
    select(paste0('fatepotential_', treatment, '_', match_name(time1), '_', match_name(time2)))
  
  fp_weighted_dist.lin <- weighted_mean_by_fatepotential(dist_mat.percell.mean, fatepotential_t1_fot_t2.lin.df)
  
  dist_df <- rbind(dist_df, data.frame(lineage = lin,
                                       dist = fp_weighted_dist.lin,
                                       n_cells.t1 = length(cell.time1),
                                       n_cells.t2 = length(cell.time2)))
}

# random sampling
dist_df.random <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(dist_df.random) <- c('lineage', 'dist.random', 'n_cells.t1', 'n_cells.t2')

for (lin in lin.common) {
  
  n.cell.time1 <- rownames(metadat.time1[metadat.time1$assigned_lineage == lin, ]) %>% length()
  n.cell.time2 <- rownames(metadat.time2[metadat.time2$assigned_lineage == lin, ]) %>% length()
  
  if(n.cell.time1 < 2 | n.cell.time2 < 2) {
    print(paste0("Not enough cells in lineage ", lin, " at time point ", time1, " or ", time2))
    next
  }
  
  mean_dist <- lapply(1: 100, function(i) {
    
    metadat.time1.low_fp <- metadat.time1[metadat.time1[[paste0('fatepotential_', treatment, '_', match_name(time1), '_', match_name(time2) )]] < 0, ]
    cell.time1.sample <- sample(rownames(metadat.time1.low_fp), n.cell.time1)
    cell.time2.sample <- sample(rownames(metadat.time2), n.cell.time2)
    
    pca.time1.lin.sample <- pca.time1[cell.time1.sample, ]
    pca.time2.lin.sample <- pca.time2[cell.time2.sample, ]
    
    
    # Calculate the pairwise distances
    dist_mat <- euclidean_dist(pca.time1.lin.sample, pca.time2.lin.sample)
    
    # Calculate row median
    dist_mat.percell.mean <- colMeans(dist_mat, na.rm = TRUE)
    
    # Get the mean distance for this lineage
    fatepotential_t1_fot_t2.lin.df <- metadat.time1[cell.time1.sample, ] %>% 
      select(paste0('fatepotential_', treatment, '_', match_name(time1), '_', match_name(time2)))
    
    fp_weighted_dist.lin <- weighted_mean_by_fatepotential(dist_mat.percell.mean, fatepotential_t1_fot_t2.lin.df)
    
    # Get the weighted distance for this lineage
    return(fp_weighted_dist.lin)
  })
  
  mean_dist.sampleMean <- mean(unlist(mean_dist))
  
  dist_df.random <- rbind(dist_df.random, data.frame(lineage = lin,
                                                     dist.random = mean_dist.sampleMean,
                                                     n_cells.t1 = n.cell.time1,
                                                     n_cells.t2 = n.cell.time2))
}


# merge dataframes
dist_df <- merge(dist_df, dist_df.random, by = c('lineage', 'n_cells.t1', 'n_cells.t2'), all = TRUE)
dist_df$adaptation.index <- dist_df$dist / dist_df$dist.random

max_val <- stats::quantile(dist_df$adaptation.index, 
                           probs = 0.99, 
                           na.rm = TRUE)
min_val <- stats::quantile(dist_df$adaptation.index,
                           probs = 0.01, 
                           na.rm = TRUE)                           
dist_df$adaptation.index_rescaled <- scales::rescale(
  pmax(pmin(dist_df$adaptation.index, 
            max_val), 
       min_val),
  to = c(0, 1)
)


write.csv(dist_df, 
          file = paste0(output_dir, 'adaptation_index_', treatment, '_', time1, '_', time2, '.csv'), 
          row.names = FALSE)


lin <- 'Lin39814'
adaptation.index <- dist_df$adaptation.index[dist_df$lineage == lin]
cells.check.time1 <- metadat.time1[metadat.time1$assigned_lineage == lin, ]
cells.check.time2 <- metadat.time2[metadat.time2$assigned_lineage == lin, ]


umap <- all_data_ft_COCL2_umap@cell.embeddings
umap <- as.data.frame(umap)
umap$cell_id <- rownames(umap)

metadat.time1$cell_id <- rownames(metadat.time1)

umap <- merge(umap, 
              metadat.time1[, c(paste0('fatepotential_', treatment, '_', match_name(time1), '_', match_name(time2)), 'cell_id')], 
              by = 'cell_id', 
              all.x = TRUE)

# umap.cells.time1 <- umap[rownames(umap) %in% rownames(cells.check.time1), ]
# umap.cells.time2 <- umap[rownames(umap) %in% rownames(cells.check.time2), ]

umap.cells.time1 <- umap[umap$cell_id %in% rownames(cells.check.time1), ]
umap.cells.time2 <- umap[umap$cell_id %in% rownames(cells.check.time2), ]

ggplot(umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  # geom_point(data = umap.cells.time1, color='red', aes(fill = 'Day10')) +
  geom_point(data = umap.cells.time1, aes(color = fatepotential_DABTRAM_d10_w5)) +
  geom_point(data = umap.cells.time2, color='darkgreen', aes(fill = 'Week5')) +
  # ggrepel::geom_text_repel(data = umap.cells.time1, aes(label = cell_id), size = 2, max.overlaps = 20) +
  # scale_fill_manual(name="Timepoints",values = c('Day10' = 'red', 'Week5' = 'blue')) +
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(lin, ' adapt. index: ', round(adaptation.index, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')

ggplot(umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  # geom_point(data = umap.cells.time1, color='red', aes(fill = 'Day10')) +
  geom_point(data = umap.cells.time1, aes(color = fatepotential_COCL2_d10_w5)) +
  geom_point(data = umap.cells.time2, color='darkgreen', aes(fill = 'Week5')) +
  # ggrepel::geom_text_repel(data = umap.cells.time1, aes(label = cell_id), size = 2, max.overlaps = 20) +
  # scale_fill_manual(name="Timepoints",values = c('Day10' = 'red', 'Week5' = 'blue')) +
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
  ggtitle(paste0(lin, ' adapt. index: ', round(adaptation.index, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')
