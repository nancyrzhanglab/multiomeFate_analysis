rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# compute nearest-neighbors
num_neigh <- 50

set.seed(10)
mat <- all_data[["pca"]]@cell.embeddings
nn_mat_rna <- RANN::nn2(mat, k = num_neigh+1)$nn.idx
nn_mat_rna <- nn_mat_rna[,-1]

set.seed(10)
mat <- all_data[["lsi"]]@cell.embeddings
nn_mat_atac <- RANN::nn2(mat, k = num_neigh+1)$nn.idx
nn_mat_atac <- nn_mat_atac[,-1]

# compute similarity
n <- ncol(all_data)
jaccard_vec <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  vec1 <- nn_mat_rna[i,]
  vec2 <- nn_mat_atac[i,]
  length(intersect(vec1, vec2))/length(unique(c(vec1, vec2)))
})

quantile(jaccard_vec)
scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(jaccard_vec^val, rank(jaccard_vec))
})
all_data$jaccard <- jaccard_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(all_data$jaccard), max(all_data$jaccard), length.out = num_color)
color_vec <- sapply(all_data$jaccard, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})
names(jaccard_vec) <- colnames(all_data)

########################################

set.seed(10)
Seurat::Idents(all_data) = "dataset"
all_data <- Seurat::RunUMAP(all_data, reduction.key="ftUMAP", reduction = "fasttopic_DABTRAM", dims=c(1:30))

# How many cells need to be seen to count as "win"?
min_cell_count <- 50
earlierday_cond_str <- "day10_DABTRAM"
later_cond <- "week5_DABTRAM"

lintab <- as.data.frame.matrix(table(all_data$assigned_lineage, all_data$dataset)) # lineage by condition table

sel.lineages <- rownames(lintab)[which(lintab[,later_cond] >= min_cell_count)]  # lineages that survived in given condition
cells.1 <- which(all_data$assigned_lineage %in% sel.lineages 
                 & all_data$assigned_posterior > 0.5
                 & all_data$dataset == earlierday_cond_str)

cells.2 <- which(!all_data$assigned_lineage %in% sel.lineages 
                 & all_data$assigned_posterior > 0.5
                 & all_data$dataset == earlierday_cond_str)

#########################################################

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_jaccard_frontrunner-", min_cell_count, ".png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = all_data@reductions[["umap"]]@cell.embeddings[,1],
     y = all_data@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Jaccard"),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "",
     cex.main = 0.5)

points(x = all_data@reductions[["umap"]]@cell.embeddings[cells.1,1],
       y = all_data@reductions[["umap"]]@cell.embeddings[cells.1,2],
       col = "white", pch = 16, cex = 3.5)
points(x = all_data@reductions[["umap"]]@cell.embeddings[cells.1,1],
       y = all_data@reductions[["umap"]]@cell.embeddings[cells.1,2],
       col = color_vec[cells.1], pch = 16, cex = 2)
graphics.off()

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_jaccard.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = all_data@reductions[["umap"]]@cell.embeddings[,1],
     y = all_data@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Synchrony"),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "",
     cex.main = 0.5)
graphics.off()


png(paste0("../../../../out/figures/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_jaccard-histogram.png"),
    height = 1500, width = 3000, units = "px", res = 500)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))

hist(all_data$jaccard[cells.1], breaks = 50, col = "gray", 
     main = paste0("Cells in ",earlierday_cond_str," whose lineage\nhas > ",min_cell_count," cells in ",later_cond),
     cex.main = 0.6,
     xlim = range(all_data$jaccard))
mean_val <- mean(all_data$jaccard[cells.1])
median_val <- median(all_data$jaccard[cells.1])
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2, lty = 1)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

hist(all_data$jaccard[cells.2], breaks = 50, col = "gray", 
     main = paste0("Cells in ",earlierday_cond_str," whose lineage\nhas < ",min_cell_count," cells in ",later_cond),
     cex.main = 0.6,
     xlim = range(all_data$jaccard))
mean_val <- mean(all_data$jaccard[cells.2])
median_val <- median(all_data$jaccard[cells.2])
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2, lty = 1)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

graphics.off()
