rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# How many cells need to be seen to count as "win"?
min_cell_count <- 50
earlierday_cond_str <- "day10_DABTRAM"
later_cond <- "week5_DABTRAM"

lintab <- as.data.frame.matrix(table(all_data$assigned_lineage, all_data$dataset)) # lineage by condition table

sel.lineages <- rownames(lintab)[which(lintab[,later_cond] >= min_cell_count)]  # lineages that survived in given condition
cells.1 <- which(all_data$assigned_lineage %in% sel.lineages 
                 & all_data$assigned_posterior > 0.5
                 & all_data$dataset == earlierday_cond_str)

# compute nearest-neighbors
set.seed(10)
mat <- all_data[["pca"]]@cell.embeddings
num_neigh <- 30
nn_mat <- RANN::nn2(mat, k = num_neigh+1)$nn.idx

# compute synchrony scores
rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

# form smoothing nearest-neighbor matrix
n <- nrow(nn_mat)
i_vec <- rep(1:n, each = ncol(nn_mat))
j_vec <- unlist(lapply(1:n, function(i){nn_mat[i,]}))
avg_mat <- Matrix::sparseMatrix(i = i_vec, 
                                j = j_vec, 
                                x = rep(1/ncol(nn_mat), length(i_vec)), 
                                dims = c(n,n))
avg_mat <- Matrix::t(avg_mat)
alignment_vec_smoothed <- as.numeric(alignment_vec %*% avg_mat)

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec_smoothed^val, rank(alignment_vec_smoothed))
})
all_data$alignment <- alignment_vec_smoothed^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(all_data$alignment), max(all_data$alignment), length.out = num_color)
color_vec <- sapply(all_data$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})
names(alignment_vec_smoothed) <- colnames(all_data)

##############

# plotting

set.seed(10)
Seurat::Idents(all_data) = "dataset"
all_data <- Seurat::RunUMAP(all_data, reduction.key="ftUMAP", reduction = "fasttopic_DABTRAM", dims=c(1:30))

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_DABTRAM_frontrunner-", min_cell_count, "_synchrony_nn.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = all_data@reductions[["umap"]]@cell.embeddings[,1],
     y = all_data@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Synchrony"),
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

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_DABTRAM_synchrony_nn.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = all_data@reductions[["umap"]]@cell.embeddings[,1],
     y = all_data@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Synchrony"),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "",
     cex.main = 0.5)
graphics.off()

# plot the histogram

cells.2 <- which(!all_data$assigned_lineage %in% sel.lineages 
                 & all_data$assigned_posterior > 0.5
                 & all_data$dataset == earlierday_cond_str)

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_DABTRAM_synchrony-histogram_nn.png"),
    height = 1500, width = 3000, units = "px", res = 500)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))

hist(all_data$alignment[cells.1], breaks = 50, col = "gray", 
     main = paste0("Cells in ",earlierday_cond_str," whose lineage\nhas > ",min_cell_count," cells in ",later_cond),
     cex.main = 0.6,
     xlim = range(all_data$alignment))
mean_val <- mean(all_data$alignment[cells.1])
median_val <- median(all_data$alignment[cells.1])
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2, lty = 1)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

hist(all_data$alignment[cells.2], breaks = 50, col = "gray", 
     main = paste0("Cells in ",earlierday_cond_str," whose lineage\nhas < ",min_cell_count," cells in ",later_cond),
     cex.main = 0.6,
     xlim = range(all_data$alignment))
mean_val <- mean(all_data$alignment[cells.2])
median_val <- median(all_data$alignment[cells.2])
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2, lty = 1)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

graphics.off()

