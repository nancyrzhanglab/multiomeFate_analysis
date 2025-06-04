rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
score_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V3/Fig6/'


remove_unassigned_cells <- TRUE

treatment <- 'DABTRAM'
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
pca <- all_data@reductions[["pca"]]@cell.embeddings

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

pca.time1 <- pca[rownames(pca) %in% rownames(metadat.time1), ]
pca.time2 <- pca[rownames(pca) %in% rownames(metadat.time2), ]

# =============================================================================
# Helper to calculate pairwise euclidean distance
# =============================================================================
euclidean_dist <- function(dimrec_t1, dimrec_t2) {
  # for each row in dimrec_t1, calculate the euclidean distance to each row in dimrec_t2
  # and return a matrix of distances
  
  # Ensure A and B are matrices
  A <- as.matrix(dimrec_t1)
  B <- as.matrix(dimrec_t2)
  
  # Compute squared row sums
  A_sq <- rowSums(A^2)
  B_sq <- rowSums(B^2)

  # print(paste0("A: ", dim(A)[1], " x ", dim(A)[2]))
  # print(paste0("B: ", dim(B)[1], " x ", dim(B)[2]))
  
  # Use matrix multiplication to compute pairwise distances
  dist_sq <- outer(A_sq, B_sq, "+") - 2 * (A %*% t(B))
  dist_sq[dist_sq < 0] <- 0  # Avoid negative due to floating point precision
  
  # Return the Euclidean distances
  sqrt(dist_sq)
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
  dist_mat.percell.mean <- rowMeans(dist_mat, na.rm = TRUE)
  
  # Get the mean distance for this lineage
  mean_dist <- mean(dist_mat.percell.mean, na.rm = TRUE)
  
  dist_df <- rbind(dist_df, data.frame(lineage = lin,
                                       dist = mean_dist,
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
  
    cell.time1.sample <- sample(rownames(metadat.time1), n.cell.time1)
    cell.time2.sample <- sample(rownames(metadat.time2), n.cell.time2)
    
    pca.time1.lin.sample <- pca.time1[cell.time1.sample, ]
    pca.time2.lin.sample <- pca.time2[cell.time2.sample, ]
    
    
    # Calculate the pairwise distances
    dist_mat <- euclidean_dist(pca.time1.lin.sample, pca.time2.lin.sample)
    
    # Calculate row median
    dist_mat.percell.mean <- rowMeans(dist_mat, na.rm = TRUE)
    
    # Get the mean distance for this lineage
    return(mean(unlist(dist_mat.percell.mean), na.rm = TRUE))
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

lin <- 'Lin129134'
adaptation.index <- dist_df$adaptation.index[dist_df$lineage == lin]
cells.check.time1 <- metadat.time1[metadat.time1$assigned_lineage == lin, ]
cells.check.time2 <- metadat.time2[metadat.time2$assigned_lineage == lin, ]

umap <- all_data_ft_COCL2_umap@cell.embeddings
umap <- as.data.frame(umap)
umap.cells.time1 <- umap[rownames(umap) %in% rownames(cells.check.time1), ]
umap.cells.time2 <- umap[rownames(umap) %in% rownames(cells.check.time2), ]

ggplot(umap, aes(x = ftCOCL2umap_1, y = ftCOCL2umap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  geom_point(data = umap.cells.time1, color='red', aes(fill = 'Day10')) +
  geom_point(data = umap.cells.time2, color='blue', aes(fill = 'Week5')) +
  scale_fill_manual(name="Timepoints",values = c('Day10' = 'red', 'Week5' = 'blue')) +
  ggtitle(paste0(lin, ' adapt. index: ', round(adaptation.index, 2))) +
  theme_classic() +
  theme(legend.position = 'bottom')

ggplot(dist_df) +
  geom_histogram(aes(x = adaptation.index), 
                 color = 'navy', fill = 'lightblue', bins = 50) +
  theme_bw()

# ggplot for histogram and desnity plot
ggplot(dist_df) +
  geom_histogram(aes(x = adaptation.index), 
                 color = 'navy', fill = 'lightblue', bins = 50) +
  geom_density(aes(x = adaptation.index), 
               color = 'navy', size = 1) +
  theme_bw() +
  xlab('Adaptation index') +
  ylab('Count') +
  ggtitle('Adaptation index distribution: CIS')


ggplot(dist_df) +
  geom_histogram(aes(x = dist),
                 color = 'navy', fill = 'blue', bins = 50, alpha = 0.5) +
  geom_histogram(aes(x = dist.random),
                 color = 'black', fill = 'gray', bins = 50, alpha = 0.5) +
  theme_bw() +
  xlab('Average Euclidean Distance') +
  ylab('Count') +
  ggtitle('Distance distribution: CIS')


p1 <- ggplot(dist_df, aes(x = log10(n_cells.t1), y = adaptation.index)) +
  geom_point() +
  stat_cor() +
  theme_bw()

p2 <- ggplot(dist_df, aes(x = log10(n_cells.t2), y = adaptation.index)) +
  geom_point() +
  stat_cor() +
  theme_bw()

ggarrange(p1, p2)

# ==============================================================================
# Compare with FP gini
# ==============================================================================
fp_gini <- read.csv(file.path(output_dir, 'fate_potential_gini_index.csv'))
fp_gini <- fp_gini[fp_gini$dataset %in% paste0(treatment, '_d10_w5'), ]

dist_df <- merge(dist_df, fp_gini, by = 'lineage', all = TRUE)
dist_df.to_plot <- dist_df %>% drop_na()

ggplot(dist_df.to_plot, aes(x = adaptation.index, y = gini_index)) +
  geom_point(aes(size = log10(n_cells.t2))) +
  geom_vline(xintercept = 0.4, linetype = 'dashed', color = 'gray') +
  # geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
  scale_color_gradient(low = 'blue', high = 'red') +
  scale_size_continuous(range = c(1, 5)) +
  xlim(0, max(dist_df.to_plot$adaptation.index, rm.na=T)) +
  ylim(0, 1) +
  theme_classic() +
  xlab('Adaptation index') +
  ylab('Fate potential gini index')

# ==============================================================================
# Compare with multi-resistant clones
# ==============================================================================
metadat.week5 <- metadat[metadat$dataset %in% c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS'), ]
lin.size.week5 <- table(metadat.week5$assigned_lineage, 
                        metadat.week5$dataset)
lin.size.week5 <- as.data.frame(lin.size.week5)
lin.size.week5 <- pivot_wider(lin.size.week5, 
                              names_from = Var2, 
                              values_from = Freq)
colnames(lin.size.week5)[1] <- 'assigned_lineage'

# sum non-zero entries per row
lin.size.week5$multi_resistant <- rowSums(lin.size.week5[, c('week5_DABTRAM', 'week5_COCL2', 'week5_CIS')] > 0)
lin.multi_resistant <- lin.size.week5[lin.size.week5$multi_resistant > 1, ] %>% pull(assigned_lineage)

dist_df.to_plot$multi_resistant_bi <- dist_df.to_plot$lineage %in% lin.multi_resistant
dist_df.to_plot <- merge(dist_df.to_plot, 
                         lin.size.week5[, c('assigned_lineage', 'multi_resistant')], 
                         by.x = 'lineage', by.y = 'assigned_lineage', all.x = TRUE)
dist_df.to_plot$multi_resistant <- as.character(dist_df.to_plot$multi_resistant)
ggplot(dist_df.to_plot, aes(x = adaptation.index, y = gini_index)) +
  geom_point(aes(size = log10(n_cells.t2), fill = multi_resistant), shape = 21) +
  # geom_vline(xintercept = quantile(dist_df.to_plot$adaptation.index, 0.5), linetype = 'dashed', color = 'gray') +
  # geom_hline(yintercept = quantile(dist_df.to_plot$gini_index, 0.5), linetype = 'dashed', color = 'gray') +
  scale_fill_manual(values = c('1' = 'gray', '2' = 'pink', '3' = 'red')) +
  scale_size_continuous(range = c(1, 15), breaks=c(0, 0.5, 1, 2, 3)) +
  ylim(0, 1) +
  theme_classic() +
  xlab('Adaptation index') +
  ylab('Fate potential gini index')

ggplot(dist_df.to_plot, aes(x = multi_resistant, y = adaptation.index)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_classic() +
  xlab('Multi-resistant clones') +
  ylab('Adaptation index') +
  ggtitle('Adaptation index distribution: CIS')
# ==============================================================================
# Compute mean fate potential
# ==============================================================================
# fp_dabtram_d10_w5 <- all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["cell_imputed_score"]]
# fp_dabtram_d10_w5 <- as.data.frame(fp_dabtram_d10_w5)
# fp_dabtram_d10_w5$cell_id <- rownames(fp_dabtram_d10_w5)
# 
# fp_dabtram_d10_w5 <- merge(fp_dabtram_d10_w5, metadat.time1[, c('cell_id', 'assigned_lineage')], by = "cell_id")
# 
# fp_dabtram_d10_w5_lineage_mean <- fp_dabtram_d10_w5 %>%
#   group_by(assigned_lineage) %>%
#   summarise(median_fp_dabtram_d10_w5 = median(fp_dabtram_d10_w5, na.rm = TRUE))
# 
# fp_lin_dabtram_d10_w5 <- all_data_fatepotential[["fatepotential_DABTRAM_d10_w5"]][["lineage_imputed_count"]]
# fp_lin_dabtram_d10_w5 <- as.data.frame(fp_lin_dabtram_d10_w5)
# fp_lin_dabtram_d10_w5$assigned_lineage <- rownames(fp_lin_dabtram_d10_w5)

# ==============================================================================
# Compare with plasticity 
# ==============================================================================
results_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task0_explore_lineage_variability_V2/'
lin_var.day10_DABTRAM <- read.csv(paste0(results_dir, 'day10_DABTRAM/lineage_variability_day10_DABTRAM_saver_sample.csv'))

dist_df.to_plot <- merge(dist_df.to_plot, 
                         lin_var.day10_DABTRAM[, c('assigned_lineage', 'normalized_avg_eud_dist_by_shuffle')], 
                         by.x = 'lineage', by.y = 'assigned_lineage', all.x = TRUE)

# ==============================================================================
# Compare with adaptaiton gene expression at day10
# ==============================================================================
scores <- read.csv(paste0(score_dir, 'UCell_scores_adaptation_gene_sets.csv'), row.names = 1)
scores$dabtram.adaptation.genes_scale_UCell <- scale(scores$dabtram.adaptation.genes_UCell)
scores$cocl2.adaptation.genes_scale_UCell <- scale(scores$cocl2.adaptation.genes_UCell)

scores.time1 <- scores[rownames(scores) %in% rownames(metadat.time1), ]
scores.time1$cell_id <- rownames(scores.time1)
scores.time2 <- scores[rownames(scores) %in% rownames(metadat.time2), ]

metadat.time1$cell_id <- rownames(metadat.time1)
metadat.time1 <- merge(metadat.time1, scores.time1, by = 'cell_id', all = TRUE)

mean.scores.time1 <- metadat.time1 %>%
  group_by(assigned_lineage) %>%
  summarise_at(vars(ends_with('UCell')), median, na.rm = TRUE)

max.scores.time1 <- metadat.time1 %>%
  group_by(assigned_lineage) %>%
  summarise_at(vars(ends_with('UCell')), max, na.rm = TRUE)

dist_df.to_plot <- merge(dist_df.to_plot, mean.scores.time1, by.x = 'lineage', by.y = 'assigned_lineage')
dist_df.to_plot <- merge(dist_df.to_plot, max.scores.time1, by.x = 'lineage', by.y = 'assigned_lineage')
# dist_df.to_plot <- merge(dist_df.to_plot, fp_dabtram_d10_w5_lineage_mean, by.x = 'lineage', by.y = 'assigned_lineage')
# dist_df.to_plot <- merge(dist_df.to_plot, fp_lin_dabtram_d10_w5, by.x = 'lineage', by.y = 'assigned_lineage')

# midpoint <- mean(dist_df.to_plot$dabtram.adaptation.genes_UCell)
midpoint <- mean(c(scores.time1$dabtram.adaptation.genes_UCell, scores.time2$dabtram.adaptation.genes_UCell))
# midpoint <- mean(scores.time2$cocl2.adaptation.genes_scale_UCell)
midpoint <- mean(dist_df.to_plot$normalized_avg_eud_dist_by_shuffle)

dist_df.to_plot <- dist_df.to_plot %>% 
  arrange(desc(fp_lin_dabtram_d10_w5))

ggplot(dist_df.to_plot, aes(x = adaptation.index , y = gini_index)) +
  geom_vline(xintercept = quantile(dist_df.to_plot$adaptation.index, 0.5), linetype = 'dashed', color = 'gray') +
  geom_hline(yintercept = quantile(dist_df.to_plot$gini_index, 0.5), linetype = 'dashed', color = 'gray') +
  geom_point(aes(size = log10(n_cells.t2), fill = dabtram.adaptation.genes_UCell), shape = 21) +
  # geom_point(aes(size = log10(n_cells.t2), fill = normalized_avg_eud_dist_by_shuffle), shape = 21) +
  # geom_point(aes(size = median_fp_dabtram_d10_w5, fill = dabtram.adaptation.genes_UCell), shape = 21) +
  # ggrepel::geom_text_repel(aes(label = lineage), size = 3, max.overlaps = 20) +
  scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
  scale_size_continuous(range = c(1, 15), breaks=c(0, 0.5, 1, 2, 3)) +
  ylim(0, 1) +
  # xlim(0, 1) +
  stat_cor() +
  theme_classic() +
  xlab('Adaptation index') +
  ylab('Fate potential gini index') +
  theme(legend.position = 'right')
ggsave(file.path(figure_dir, 'adaptation_index_vs_gini_index_DABTRAM_d10_w5.pdf'),
       width = 7.5, height = 5)

ggplot(dist_df.to_plot, aes(x = adaptation.index , y = gini_index)) +
  geom_vline(xintercept = 0.5, linetype = 'dashed', color = 'gray') +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'gray') +
  geom_point(aes(size = log10(n_cells.t2), fill = cocl2.adaptation.genes_UCell), shape = 21) +
  # geom_point(aes(size = median_fp_dabtram_d10_w5, fill = dabtram.adaptation.genes_UCell), shape = 21) +
  # ggrepel::geom_text_repel(aes(label = lineage), size = 3, max.overlaps = 20) +
  scale_fill_gradient2(low = '#309898', mid = 'bisque', high = '#CB0404', midpoint = midpoint) +
  scale_size_continuous(range = c(1, 15), breaks=c(0, 0.5, 1, 2, 3)) +
  ylim(0.45, 1) +
  xlim(0.5, 1.2) +
  theme_classic() +
  xlab('Adaptation index') +
  ylab('Fate potential gini index') +
  theme(legend.position = 'right')


lin <- 'Lin130951'
cells.check.time1 <- metadat.time1[metadat.time1$assigned_lineage == lin, ]
cells.check.time2 <- metadat.time2[metadat.time2$assigned_lineage == lin, ]

umap <- all_data_ft_DABTRAM_umap@cell.embeddings
umap <- as.data.frame(umap)
umap.cells.time1 <- umap[rownames(umap) %in% cells.check.time1$cell_id, ]
umap.cells.time2 <- umap[rownames(umap) %in% rownames(cells.check.time2), ]

umap.cells.time1 <- merge(umap.cells.time1, scores.time1, by.x = 'row.names', by.y = 'cell_id')
umap.cells.time1 <- merge(umap.cells.time1, fp_dabtram_d10_w5, by.x = 'Row.names', by.y = 'cell_id')

scores.time1.all <- merge(scores.time1, umap.cells.time1[, c('Row.names', 'ftDABTRAMumap_1', 'ftDABTRAMumap_2')], by.x = 'cell_id', by.y = 'Row.names', all = T)

ggplot(umap, aes(x = ftDABTRAMumap_1, y = ftDABTRAMumap_2)) +
  geom_point(size = 0.5, color = '#D3D3D3') +
  geom_point(data = umap.cells.time1, aes(color = dabtram.adaptation.genes_UCell)) +
  geom_point(data = scores.time1.all, aes(color = dabtram.adaptation.genes_UCell), size = 0, alpha = 0) +
  # geom_point(data = umap.cells.time2, color='black') +
  scale_color_gradient(low = 'blue', high = 'red') +
  theme_classic()
