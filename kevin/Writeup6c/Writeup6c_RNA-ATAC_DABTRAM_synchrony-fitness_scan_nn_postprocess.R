rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan_nn.RData")

rsquare_value_fit <- sapply(1:ncol(pred_diff_mat), function(j){
  tmp_df <- data.frame(y = pred_diff_mat[,j], x = log1p(future_num_vec_smoothed[order_vec]))
  
  tmp_lm <- stats::lm(y ~ x, data = tmp_df)
  summary(tmp_lm)$r.squared
})
slope_fit <- sapply(1:ncol(pred_diff_mat), function(j){
  tmp_df <- data.frame(y = pred_diff_mat[,j], x = log1p(future_num_vec_smoothed[order_vec]))
  
  tmp_lm <- stats::lm(y ~ x, data = tmp_df)
  stats::coef(tmp_lm)["x"]
})

round(quantile(rsquare_value_fit),2)
round(quantile(slope_fit),2)

gene_idx <- intersect(which(rsquare_value_fit >= 0.9),
                      which(slope_fit >= 0.25))
length(gene_idx)
