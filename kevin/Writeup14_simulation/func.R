.preprocess_rna <- function(rna_mat, treatment){
  cell_idx <- which(all_data$dataset == treatment)
  rna_mat <- rna_mat[cell_idx,]
  rna_mat <- pmin(rna_mat, 10)
  
  mean_vec <- colMeans(rna_mat)
  sd_vec <- apply(rna_mat, 2, stats::sd)
  rm_idx <- which(sd_vec <= 1e-3)
  if(length(rm_idx) > 0){
    rna_mat <- rna_mat[,-rm_idx,drop = FALSE]
  }
  cell_features <- scale(rna_mat)
  cell_features
}

.compute_mean_total_cells <- function(cell_features,
                                      coefficient_intercept,
                                      coefficient_vec){
  sum(exp((cell_features %*% coefficient_vec) + coefficient_intercept))
}

##################

.plot_mean_variance <- function(filename,
                                simulation_res){
  
  median_vec <- simulation_res$summary_mat["median",]
  range_vec <- simulation_res$summary_mat["range",]
  future_vec <- simulation_res$summary_mat["future_size",]
  
  grDevices::png(filename = filename,
                 width = 15, 
                 height = 5, 
                 units = "in", 
                 res = 300)
  par(mfrow = c(1,3))
  
  graphics::plot(median_vec, 
                 future_vec, 
                 pch = 16,
                 xlab = "Median (per lineage)",
                 ylab = "Future size (per lineage)",
                 main = paste0("Future size vs. current median potential\n",
                               "Corr: ", round(cor(median_vec, future_vec), 2)))
  graphics::plot(range_vec,
                 future_vec,
                 pch = 16,
                 xlab = "Range (per lineage)",
                 ylab = "Future size (per lineage)",
                 main = paste0("Future size vs. current range potential\n",
                               "Corr: ", round(cor(range_vec, future_vec), 2)))
  
  graphics::plot(x = median_vec,
                 y = range_vec, 
                 main = paste0("Corr: ", round(cor(range_vec, median_vec), 2)),
                 xlab = "Median (per lineage)", 
                 ylab = "Range (per lineage)",
                 pch = 16, 
                 asp = TRUE)
  
  grDevices::graphics.off()
}