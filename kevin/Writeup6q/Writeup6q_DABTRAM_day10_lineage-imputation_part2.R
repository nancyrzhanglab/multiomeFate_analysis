rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

treatment <- "DABTRAM"
day_early <- "day10"
day_later <- "week5"
day_early_full <- paste0(day_early, "_", treatment)
day_later_full <- paste0(day_later, "_", treatment)

load(paste0("../../../../out/kevin/Writeup6q/Writeup6q_", treatment, "_", day_early, "_lineage-imputation.RData"))

###################

lambda_sequence <- loocv_fit_list[[1]]$lambda_sequence
stopifnot(sum(abs(lambda_sequence - loocv_fit_list[[2]]$lambda_sequence)) <= 1e-4)

# for each lineage being left out, compute the sequence of fits
loocv_lineages <- names(loocv_fit_list)
loocv_len <- length(loocv_lineages)
loocv_idx_list <- lapply(loocv_lineages, function(lineage){
  which(cell_lineage == lineage)
})
names(loocv_idx_list) <- loocv_lineages

# just run a quick check
for(lineage_name in loocv_lineages){
  stopifnot(length(loocv_idx_list[[lineage_name]]) == lineage_current_count[lineage_name])
  stopifnot(lineage_future_count[lineage_name] >= 20)
}

res_mat <- matrix(NA, nrow = loocv_len, ncol = length(lambda_sequence))
rownames(res_mat) <- loocv_lineages
colnames(res_mat) <- lambda_sequence

for(lineage in loocv_lineages){
  print(paste0("Working on lineage ", lineage))
  for(j in 1:length(lambda_sequence)){
    lambda <- lambda_sequence[j]
    
    cell_features_test <- cell_features[loocv_idx_list[[lineage]],,drop = F]
    cell_lineage_test <- cell_lineage[loocv_idx_list[[lineage]]]
    lineage_future_count_test <- lineage_future_count[lineage]
    
    res_mat[lineage,j] <- multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_test,
                                                                cell_lineage = cell_lineage_test,
                                                                coefficient_vec = loocv_fit_list[[lineage]]$fit_list[[j]]$coefficient_vec,
                                                                lineage_future_count = lineage_future_count_test,
                                                                lambda = lambda)
  }
}

loglik_mean <- Matrix::colMeans(res_mat)
