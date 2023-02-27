rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan.RData")

order_vec <- order(future_num_vec, decreasing = F)
rsquare_value_fit <- sapply(1:ncol(pred_diff_mat), function(j){
  tmp_df <- data.frame(y = pred_diff_mat[,j], x = log1p(future_num_vec[order_vec]))
  
  tmp_lm <- stats::lm(y ~ x, data = tmp_df)
  summary(tmp_lm)$r.squared
})
slope_fit <- sapply(1:ncol(pred_diff_mat), function(j){
  tmp_df <- data.frame(y = pred_diff_mat[,j], x = log1p(future_num_vec[order_vec]))
  
  tmp_lm <- stats::lm(y ~ x, data = tmp_df)
  stats::coef(tmp_lm)["x"]
})

round(quantile(rsquare_value_fit),2)
round(quantile(slope_fit),2)

gene_idx <- intersect(which(rsquare_value_fit >= 0.9),
                      which(slope_fit >= 0.25))
length(gene_idx)

###################

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))
diff_mat <- (rna_mat2 - atac_mat2)/rna_mat2
colnames(diff_mat) <- colnames(rna_mat2)
diff_mat <- abs(diff_mat)
diff_mat <- pmin(diff_mat, 15)

###################

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan_genes.png"),
    height = 1500, width = 3000, units = "px", res = 300)
n <- nrow(pred_diff_mat)
par(mar = c(4,4,4,0.5), mfrow = c(3,5))

for(i in gene_idx[1:30]){
  plot(x = log1p(future_num_vec[order_vec]),
       y = diff_mat[,i], pch = 16, col = rgb(0.75, 0.75, 0.75, 0.15),
       xlab = "Future fitness", ylab = "Synchrony", 
       main = colnames(pred_diff_mat)[i],
       ylim = c(min(diff_mat[,i]), 1.5*max(pred_diff_mat[,i])))
  lines(x = log1p(future_num_vec[order_vec]),
        y = pred_diff_mat[,i], col = "white", lwd = 6)
  lines(x = log1p(future_num_vec[order_vec]),
        y = pred_diff_mat[,i], col = "black", lwd = 4)
}
graphics.off()