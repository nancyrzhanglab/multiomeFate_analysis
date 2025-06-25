rm(list=ls())
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)

set.seed(123)

data_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/FINAL/'
score_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task5_FP_for_d10_adapting_cells/'
output_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task9_quantify_adaptation/'
figure_dir <- '/Users/emiliac/Dropbox/Thesis/Lineage_trace/figures/MultiomeFate_V5/Figure5/'


remove_unassigned_cells <- TRUE

treatment <- 'COCL2'
time1 <- 'day10'
time2 <- 'week5'

theme_Publication<- function(base_size=12, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.box.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "mm"),
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
            strip.background=element_rect(colour="#F0F0F0",fill="#F0F0F0"),
            strip.text = element_text(face="bold")
    ))
}

dataset_colors <- c(day0 = "gray",
                    day10_CIS = "#FBD08C",
                    day10_COCL2 = "#6DC49C",
                    day10_DABTRAM = "#9D85BE",
                    week5_CIS = "#C96D29",
                    week5_COCL2 = "#0F8241",
                    week5_DABTRAM = "#623594")
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
colnames(dist_df) <- c('cell_id', 'dist', 'n_cells.t1', 'n_cells.t2')

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
  
  
  dist_df <- rbind(dist_df, data.frame(cell_id = cell.time1,
                                       dist = dist_mat.percell.mean))
}

# ==============================================================================
# Compare fate potential
# ==============================================================================
fp_d10_w5 <- all_data_fatepotential[[paste0("fatepotential_", treatment, "_d10_w5")]][["cell_imputed_score"]]
fp_d10_w5 <- as.data.frame(fp_d10_w5)
fp_d10_w5$cell_id <- rownames(fp_d10_w5)

comp_df <- merge(dist_df, fp_d10_w5, by = 'cell_id')

ggplot(comp_df, aes(x = dist, y = fp_d10_w5)) +
  geom_point(size = 0.5, color = dataset_colors[paste0(time1, '_', treatment)]) +
  stat_cor() +
  labs(title = paste0("Distance between ", time1, " and ", time2, " within lineages"),
       x = "Euclidean distance",
       y = "Growth potential (day 10 to week 5)") +
  theme_Publication()

ggsave(paste0(figure_dir, 'Supp_adaptation_distance_vs_growth_potential_', treatment, '_', time1, '_', time2, '.pdf'),
       width = 5, height = 5)


