rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup4f/Writeup4f_barcode_estimation.RData")
source("barcode_estimation.R")

quantile(barcode_res$beta0)
quantile(barcode_res$beta1)

png("../../../../out/figures/Writeup4f/Writeup4f_estimated_betas.png",
    height = 1500, width = 3000, units = "px", res = 300)
par(mfrow = c(1,2))
plot(log10(barcode_res$beta0), 
     log10(barcode_res$beta1), 
     xlab = "Beta 0 (Log10)", ylab = "Beta 1 (Log10)",
     pch = 16, col = rgb(0.5, 0.5, 0.5, 0.2),
     main = "Estimated betas, one per lineage", asp = T)
hist(log10(barcode_res$beta1), xlab = "Beta 1 (Log10)", main = "")
graphics.off()


n <- ncol(lin_mat); p <- nrow(lin_mat)
num_high_posterior_per_cell <- sapply(1:n, function(i){
  length(which(barcode_res$posterior_mat1[,i] > 0.5))
})
table(num_high_posterior_per_cell)
idx <- which.max(num_high_posterior_per_cell)
idx
round(barcode_res$posterior_mat1[which(barcode_res$posterior_mat1[,idx] > 0.5),idx],3)
lin_mat[which(barcode_res$posterior_mat1[,idx] > 0.5),idx]

nonzero_list <- sapply(1:n, function(i){
  .nonzero_col(lin_mat, col_idx = i, bool_value = F)
})

cell_assignments <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  vec <- barcode_res$posterior_mat1[,i]
  idx <- which(vec >= 0.5)
  if(length(idx) == 0) return(numeric(0))
  max_idx <- which.max(vec)
  if(length(which(vec == vec[max_idx])) > 1) return(numeric(0))
  if(max_idx %in% nonzero_list[[i]]) {
    return(c(max_idx, i))
  } else {
    return(numeric(0))
  }
})
dominant_assignment_pairs <- do.call(rbind, cell_assignments)
dominant_mat <- Matrix::sparseMatrix(i = dominant_assignment_pairs[,1],
                                     j = dominant_assignment_pairs[,2],
                                     x = rep(1, nrow(dominant_assignment_pairs)),
                                     dims = c(p,n))
length(which(Matrix::colSums(dominant_mat) == 0))
idx <- which(Matrix::colSums(dominant_mat) == 0)
table(dataset_vec[idx])

dominant_assignment_pairs <- lapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  idx <- .nonzero_col(lin_mat, col_idx = i, bool_value = F)
  val <- .nonzero_col(lin_mat, col_idx = i, bool_value = T)
  
  if(length(idx) == 0) return(numeric(0))
  if(length(idx) == 1) return(c(idx[1],i))
  
  val_sorted <- sort(val, decreasing = T)
  if(val_sorted[1] <= 2*val_sorted[2]) return(numeric(0))
  return(c(idx[which.max(val)],i))
})
dominant_assignment_pairs <- do.call(rbind, dominant_assignment_pairs)
dominant_mat_naive <- Matrix::sparseMatrix(i = dominant_assignment_pairs[,1],
                                           j = dominant_assignment_pairs[,2],
                                           x = rep(1, nrow(dominant_assignment_pairs)),
                                           dims = dim(lin_mat))
length(which(Matrix::colSums(dominant_mat_naive) == 0))
idx <- which(Matrix::colSums(dominant_mat_naive) == 0)
table(dataset_vec[idx])

###################################

cell_assignments <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  vec <- barcode_res$posterior_mat1[,i]
  idx <- which(vec > 0.5)
  if(length(idx) != 1) return(numeric(0))
  max_idx <- which.max(vec)
  return(c(max_idx, i))
})
dominant_assignment_pairs <- do.call(rbind, cell_assignments)
dominant_mat_stringent <- Matrix::sparseMatrix(i = dominant_assignment_pairs[,1],
                                     j = dominant_assignment_pairs[,2],
                                     x = rep(1, nrow(dominant_assignment_pairs)),
                                     dims = c(p,n))
cell_idx <- which(Matrix::colSums(dominant_mat_stringent) > 0)
cell_idx2 <- setdiff(cell_idx, which(Matrix::colSums(dominant_mat_naive) > 0))
table(dataset_vec[cell_idx2])
naive_idx <- cell_idx2[which(dataset_vec[cell_idx2] == "day0")]
max_val_vec <- sapply(naive_idx, function(i){
  val <- .nonzero_col(lin_mat, col_idx = i, bool_value = T)
  max(val)
})
table(max_val_vec)