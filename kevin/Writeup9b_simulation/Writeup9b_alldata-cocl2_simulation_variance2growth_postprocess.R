rm(list=ls())
load("~/project/Multiome_fate/out/kevin/Writeup9b/Writeup9b_day10-COCL2_simulation_variance2growth.RData")

coefficient_intercept <- simulation_res_full$coefficient_intercept
embedding_coefficient_vec <- simulation_res_full$embedding_coefficient_vec

est_coefficient_vec <- final_fit_full$coefficient_vec
est_coefficient_intercept <- est_coefficient_vec["Intercept"]
est_coefficient_vec <- est_coefficient_vec[setdiff(names(est_coefficient_vec), "Intercept")]

####

lineage_future_size <- simulation_res_full$lineage_future_size
cell_fate_potential_truth <- simulation_res_full$cell_fate_potential_truth

lineage_imputed_count <- final_fit_full$lineage_imputed_count[names(lineage_future_size)]
cell_imputed_score <- final_fit_full$cell_imputed_score
if(length(names(cell_imputed_score)) == 0){
  names(cell_imputed_score) <- names(cell_fate_potential_truth)
}
cell_imputed_score <- cell_imputed_score[names(cell_fate_potential_truth)]

png(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9b/Writeup9b_alldata-cocl2_simulation_variance2growth_results.png"),
    width = 15, height = 5, units = "in", res = 300)
par(mfrow = c(1,3))
plot(x = embedding_coefficient_vec,
     y = est_coefficient_vec,
     xlab = "True coefficients (Full)",
     ylab = "Estimated coefficients",
     asp = TRUE, pch = 16)

plot(x = lineage_future_size,
     y = lineage_imputed_count,
     asp = TRUE, pch = 16,
     xlab = "Observed lineage size (Full)",
     ylab = "Fitted lineage size",
     main = paste0("Correlation: ", round(stats::cor(
       lineage_future_size,
       lineage_imputed_count
     ), 2))
)

plot(x = cell_fate_potential_truth,
     y = cell_imputed_score,
     asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.3),
     xlab = "True cell potential (Full)",
     ylab = "Estimated cell potential",
     main = paste0("Correlation: ", round(stats::cor(
       cell_fate_potential_truth,
       cell_imputed_score
     ), 2))
)

graphics.off()

######################

coefficient_intercept <- simulation_res_pca$coefficient_intercept
embedding_coefficient_vec <- simulation_res_pca$embedding_coefficient_vec

est_coefficient_vec <- final_fit_pca$coefficient_vec
est_coefficient_intercept <- est_coefficient_vec["Intercept"]
est_coefficient_vec <- est_coefficient_vec[setdiff(names(est_coefficient_vec), "Intercept")]

####

lineage_future_size <- simulation_res_pca$lineage_future_size
cell_fate_potential_truth <- simulation_res_pca$cell_fate_potential_truth

lineage_imputed_count <- final_fit_pca$lineage_imputed_count[names(lineage_future_size)]
cell_imputed_score <- final_fit_pca$cell_imputed_score
if(length(names(cell_imputed_score)) == 0){
  names(cell_imputed_score) <- names(cell_fate_potential_truth)
}
cell_imputed_score <- cell_imputed_score[names(cell_fate_potential_truth)]

png(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9b/Writeup9b_alldata-cocl2_pca_simulation_variance2growth_results.png"),
    width = 15, height = 5, units = "in", res = 300)
par(mfrow = c(1,3))
plot(x = embedding_coefficient_vec,
     y = est_coefficient_vec,
     xlab = "True coefficients (PCA)",
     ylab = "Estimated coefficients",
     asp = TRUE, pch = 16)

plot(x = lineage_future_size,
     y = lineage_imputed_count,
     asp = TRUE, pch = 16,
     xlab = "Observed lineage size (PCA)",
     ylab = "Fitted lineage size",
     main = paste0("Correlation: ", round(stats::cor(
       lineage_future_size,
       lineage_imputed_count
     ), 2))
)

plot(x = cell_fate_potential_truth,
     y = cell_imputed_score,
     asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.3),
     xlab = "True cell potential (PCA)",
     ylab = "Estimated cell potential",
     main = paste0("Correlation: ", round(stats::cor(
       cell_fate_potential_truth,
       cell_imputed_score
     ), 2))
)

graphics.off()


