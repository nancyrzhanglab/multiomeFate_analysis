rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan_nn.RData")

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))
diff_mat <- (rna_mat2 - atac_mat2)/rna_mat2
colnames(diff_mat) <- colnames(rna_mat2)
diff_mat <- abs(diff_mat)
diff_mat <- pmin(diff_mat, 15)

x_vec <- log1p(future_num_vec_smoothed[order_vec])
names(coef_vec) <- colnames(diff_mat)
gene_idx <- order(coef_vec, decreasing = T)[1:15]
names(coef_vec)[gene_idx]

pred_diff_mat <- sapply(gene_idx, function(j){
  print(j)
  
  tmp_df <- data.frame(y = diff_mat[,j], x = x_vec)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_diff_mat) <- names(coef_vec)[gene_idx]

###################

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan_nn_genes.png"),
    height = 1500, width = 3000, units = "px", res = 300)
n <- nrow(pred_diff_mat)
par(mar = c(4,4,4,0.5), mfrow = c(3,5))

for(i in 1:ncol(pred_diff_mat)){
  plot(x = x_vec,
       y = diff_mat[,gene_idx[i]], pch = 16, col = rgb(0.75, 0.75, 0.75, 0.15),
       xlab = "Future fitness", ylab = "Synchrony", 
       main = colnames(pred_diff_mat)[i])
  lines(x = x_vec,
        y = pred_diff_mat[,i], col = "white", lwd = 6)
  lines(x = x_vec,
        y = pred_diff_mat[,i], col = "black", lwd = 4)
}
graphics.off()

