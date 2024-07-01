# this simulation is for similar means but the different ranges

rm(list=ls())
library(Seurat)
library(multiomeFate)

load("../../out/Writeup9b/Writeup9b_simulation_day10-COCL2.RData")

set.seed(10)
rna_mat <- SeuratObject::LayerData(all_data, 
                                   layer = "data", 
                                   assay = "RNA",
                                   features = Seurat::VariableFeatures(all_data))
rna_mat <- Matrix::t(rna_mat)

early_idx <- which(all_data$dataset == "day10_COCL2")
rna_mat <- rna_mat[early_idx,]

gene_sparsity_vec <- sapply(1:ncol(rna_mat), function(j){
  length(multiomeFate:::.nonzero_col(rna_mat,
                                     col_idx = j,
                                     bool_value = FALSE))
})
names(gene_sparsity_vec) <- colnames(rna_mat)
gene_ordering <- names(gene_sparsity_vec)[order(gene_sparsity_vec, decreasing = TRUE)]
rna_mat <- as.matrix(rna_mat)
mean_vec <- colMeans(rna_mat)
sd_vec <- apply(rna_mat, 2, stats::sd)
rm_idx <- which(sd_vec <= 1e-3)
if(length(rm_idx) > 0){
  rna_mat <- rna_mat[,-rm_idx,drop = FALSE]
}
cell_features <- scale(rna_mat)
svd_res <- irlba::irlba(cell_features, nv = 30)
svd_features <- sweep(svd_res$u, MARGIN = 2, STATS = svd_res$d, FUN = '*')
rownames(svd_features) <- rownames(rna_mat)
colnames(svd_features) <- paste0("svd_", 1: ncol(svd_features))

coefficient_intercept <- 0.7
svd_coefficient_vec <- rep(0, ncol(svd_features))
names(svd_coefficient_vec) <- colnames(svd_features)
svd_coefficient_vec[1:5] <- seq(0.1, 0.04, length.out = 5)
coefficient_vec <- as.numeric(svd_res$v %*% svd_coefficient_vec)
names(coefficient_vec) <- colnames(cell_features)

# double-check the fate potentials are ok
tmp2 <- exp((svd_features %*% svd_coefficient_vec) + coefficient_intercept)
sum(tmp2)

original_tmp <- exp((cell_features %*% coefficient_vec) + coefficient_intercept)
sum(original_tmp)
# plot(tmp, tmp2, main = paste0("Corr: ", round(cor(tmp, tmp2), 2)))

num_lineages <- 50

############################
# not much to change after this line

set.seed(10)
simulation_res_pca <- multiomeFate:::generate_simulation_plastic(
  embedding_mat = svd_features,
  bool_add_randomness = TRUE,
  coefficient_intercept = coefficient_intercept,
  embedding_coefficient_vec = svd_coefficient_vec,
  num_lineages = num_lineages,
  lineage_mean_spread = 1, 
  lineage_sd_spread = NA,
  verbose = 3
)

par(mfrow = c(1,3))
median_vec <- simulation_res_pca$summary_mat["median",]
range_vec <- simulation_res_pca$summary_mat["range",]
future_vec <- simulation_res_pca$summary_mat["future_size",]
plot(median_vec, 
     future_vec, 
     pch = 16,
     xlab = "Median (per lineage)",
     ylab = "Future size (per lineage)",
     main = paste0("Future size vs. current median potential\n",
                   "Corr: ", round(cor(median_vec, future_vec), 2)))
plot(range_vec,
     future_vec,
     pch = 16,
     xlab = "Range (per lineage)",
     ylab = "Future size (per lineage)",
     main = paste0("Future size vs. current range potential\n",
                   "Corr: ", round(cor(range_vec, future_vec), 2)))
plot(x = median_vec,
     y = range_vec, 
     main = paste0("Corr: ", round(cor(range_vec, median_vec), 2)),
     xlab = "Median (per lineage)", ylab = "Range (per lineage)",
     pch = 16, asp = TRUE)

