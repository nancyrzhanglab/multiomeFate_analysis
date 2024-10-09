.form_simulation_seurat_fate <- function(final_fit,
                                         simulation_res){
  stopifnot(length(names(simulation_res$lineage_assignment)) > 0,
            length(names(final_fit$cell_imputed_score)) > 0,
            all(names(final_fit$cell_imputed_score) == names(simulation_res$lineage_assignment)))
  
  # compute the cell fates
  cell_features <- cbind(1, simulation_res$embedding_mat)
  true_coefficient_vec <- c(simulation_res$coefficient_intercept, simulation_res$coefficient_vec)
  names(true_coefficient_vec) <- c("Intercept", colnames(simulation_res$embedding_mat))
  fatepotential_true <- as.numeric(cell_features %*% true_coefficient_vec)
  names(fatepotential_true) <- rownames(cell_features)
  
  metadata <- data.frame(assigned_lineage = simulation_res$lineage_assignment,
                         fatepotential = final_fit$cell_imputed_score,
                         fatepotential_true = fatepotential_true)
  rownames(metadata) <- rownames(simulation_res$embedding_mat)
  
  # create the seurat object
  count_mat <- Matrix::Matrix(0, 
                              ncol = nrow(simulation_res$embedding_mat),
                              nrow = 10)
  colnames(count_mat) <- rownames(simulation_res$embedding_mat)
  rownames(count_mat) <- paste0("Gene", 1:10)
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = count_mat, 
                                           meta.data = metadata)
  
  seurat_obj[["fasttopic"]] <- Seurat::CreateDimReducObject(simulation_res$embedding_mat)
  
  # put in the  lineage sizes
  seurat_obj@misc$lineage_imputed_count <- final_fit$lineage_imputed_count
  seurat_obj@misc$lineage_observed_count <- simulation_res$lineage_future_size
  
  # put in the coefficients
  seurat_obj@misc$true_coefficient_vec <- true_coefficient_vec
  seurat_obj@misc$est_coefficient_vec <- final_fit$coefficient_vec
  
  seurat_obj
}

##################

.plot_cellFateScatterplot <- function(fatepotential_true,
                                      fatepotential,
                                      title = ""){
  
  stopifnot(length(names(fatepotential_true)) > 0,
            length(names(fatepotential)) > 0,
            all(names(fatepotential) == names(fatepotential_true)))
  
  n <- length(fatepotential)
  df <- data.frame(fatepotential = fatepotential,
                   fatepotential_true = fatepotential_true,
                   name = names(fatepotential))
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = fatepotential_true, 
                                            y = fatepotential))
  plot1 <- plot1 + ggplot2::geom_point()
  plot1 <- plot1 + ggplot2::ggtitle(paste0(
    title,
    "\nCorr:", round(stats::cor(fatepotential, fatepotential_true), 2))
  )
  plot1 <- plot1 + ggplot2::xlab("True fate potential (Log10)") + ggplot2::ylab("Predicted fate potential (Log10)")
  plot1 <- plot1 + Seurat::NoLegend() 
  
  plot1
}