rm(list=ls())

load("~/project/Multiome_fate/out/kevin/Writeup9b/Writeup9b_day10-COCL2_simulation_variance2growth.RData")

embedding_coefficient_vec <- simulation_res$embedding_coefficient_vec
fitted_coefficient_vec <- final_fit$coefficient_vec[names(embedding_coefficient_vec)]

lineage_imputed_count <- final_fit$lineage_imputed_count
cell_imputed_score <- final_fit$cell_imputed_score

png(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9b/Writeup9b_alldata-cocl2_simulation_variance2growth_results.png"),
    width = 15, height = 5, units = "in", res = 300)
par(mfrow = c(1,3))
plot(x = embedding_coefficient_vec,
     y = fitted_coefficient_vec,
     xlab = "True coefficients",
     ylab = "Estimated coefficients",
     asp = TRUE, pch = 16)

plot(x = simulation_res$lineage_future_size,
     y = lineage_imputed_count[names(simulation_res$lineage_future_size)],
     asp = TRUE, pch = 16,
     xlab = "Observed lineage size",
     ylab = "Fitted lineage size",
     main = paste0("Correlation: ", round(stats::cor(
       simulation_res$lineage_future_size,
       lineage_imputed_count[names(simulation_res$lineage_future_size)]
     ), 2))
)

plot(x = simulation_res$cell_fate_potential_truth,
     y = cell_imputed_score[names(simulation_res$cell_fate_potential_truth)],
     asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.3),
     xlab = "True cell potential",
     ylab = "Estimated cell potential",
     main = paste0("Correlation: ", round(stats::cor(
       simulation_res$cell_fate_potential_truth,
       cell_imputed_score[names(simulation_res$cell_fate_potential_truth)]
     ), 2))
)
graphics.off()
