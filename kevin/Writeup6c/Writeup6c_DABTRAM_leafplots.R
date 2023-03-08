rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_idx <- intersect(intersect(which(all_data$dataset == "day10_DABTRAM"),
                                which(all_data$assigned_posterior > 0.5)),
                      which(!is.na(all_data$assigned_lineage)))


lineage_vec <- all_data$assigned_lineage[cell_idx]
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
future_num_vec <- sapply(lineage_vec, function(lineage){
  tab_mat[lineage,"week5_DABTRAM"]
})
future_num_vec <- log1p(future_num_vec)

num_colors <- 20
color_palette <- grDevices::colorRampPalette(c("gray", "coral2"))(num_colors)
color_breaks <- seq(min(future_num_vec), max(future_num_vec), length.out = num_colors)
color_vec <- sapply(future_num_vec, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

########

rna_common <- multiSVD_obj$common_mat_1
rna_all <- multiSVD_obj$common_mat_1 + multiSVD_obj$distinct_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

p <- ncol(atac_pred)
cor_vec <- sapply(1:p, function(j){
  stats::cor(rna_all[cell_idx,j], atac_pred[cell_idx,j])
})
quantile(cor_vec)
gene_idx <- order(cor_vec, decreasing = T)[1:25]

set.seed(10)
rng_idx <- sample(1:length(cell_idx))

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_DABTRAM_leafplot_genes_highest-correlation.png"),
    height = 3000, width = 3000, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(5,5))
for(j in gene_idx){
  plot(x = atac_pred[cell_idx[rng_idx],j],
       y = rna_all[cell_idx[rng_idx],j], pch = 16, col = color_vec[rng_idx],
       xlab = "Mapped ATAC value", ylab = "RNA value", 
       main = paste0(colnames(rna_all)[j], "\nCorr: ", round(cor_vec[j],2)),
       cex.main = 0.7)
}
graphics.off()

################

rna_common <- multiSVD_obj$common_mat_1[cell_idx,]
svd_1 <- irlba::irlba(rna_common, nv = 49)
atac_common <- multiSVD_obj$common_dimred_2[cell_idx,]
svd_2 <- svd(atac_common)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

p <- ncol(atac_pred)
cor_vec <- sapply(1:p, function(j){
  stats::cor(rna_common[,j], atac_pred[,j])
})
quantile(cor_vec)
gene_idx <- order(cor_vec, decreasing = T)[1:25]


