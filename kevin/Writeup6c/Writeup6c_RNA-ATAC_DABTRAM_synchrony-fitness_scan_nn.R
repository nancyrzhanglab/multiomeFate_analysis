rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)


######################

cell_idx <- intersect(intersect(which(all_data$dataset == "day10_DABTRAM"),
                                which(all_data$assigned_posterior > 0.5)),
                      which(!is.na(all_data$assigned_lineage)))
lineage_vec <- all_data$assigned_lineage[cell_idx]

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
future_num_vec <- sapply(lineage_vec, function(lineage){
  tab_mat[lineage,"week5_DABTRAM"]
})

# compute nearest-neighbors, just among the valid day10 cells
set.seed(10)
mat <- all_data[["pca"]]@cell.embeddings[cell_idx,]
num_neigh <- 30
nn_mat <- RANN::nn2(mat, k = num_neigh+1)$nn.idx

# form smoothing nearest-neighbor matrix
n <- nrow(nn_mat)
i_vec <- rep(1:n, each = ncol(nn_mat))
j_vec <- unlist(lapply(1:n, function(i){nn_mat[i,]}))
avg_mat <- Matrix::sparseMatrix(i = i_vec, 
                                j = j_vec, 
                                x = rep(1/ncol(nn_mat), length(i_vec)), 
                                dims = c(n,n))
avg_mat <- Matrix::t(avg_mat)
future_num_vec_smoothed <- as.numeric(future_num_vec %*% avg_mat)

# compute alignment
rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

# for a gene, compute the regression line regressing alignment onto future fitness across all the cells
target_genes <- colnames(rna_common)

order_vec <- order(future_num_vec_smoothed, decreasing = F)
rna_mat <- rna_common[cell_idx[order_vec],target_genes]
atac_mat <- atac_pred[cell_idx[order_vec],target_genes]

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))
diff_mat <- (rna_mat2 - atac_mat2)/rna_mat2
colnames(diff_mat) <- colnames(rna_mat2)
diff_mat <- abs(diff_mat[,target_genes])
diff_mat <- pmin(diff_mat, 15)

x_vec <- log1p(future_num_vec_smoothed[order_vec])

p <- ncol(diff_mat); n <- nrow(diff_mat)
pred_diff_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = diff_mat[,j], x = x_vec)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_diff_mat) <- colnames(diff_mat)

coef_vec <- sapply(1:ncol(pred_diff_mat), function(j){
  tmp_df <- data.frame(y = diff_mat[,j], x = x_vec)
  
  tmp_lm <- stats::lm(y ~ x, data = tmp_df)
  # summary(tmp_lm)$r.squared
  stats::coef(tmp_lm)["x"]
})

save(pred_diff_mat, coef_vec, nn_mat, order_vec,
     rna_mat, atac_mat, future_num_vec, future_num_vec_smoothed,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan_nn.RData")

