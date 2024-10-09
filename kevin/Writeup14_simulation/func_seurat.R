.form_simulation_seurat_fate <- function(final_fit,
                                         simulation_res){
  stopifnot(length(names(simulation_res$lineage_assignment)) > 0,
            length(names(final_fit$cell_imputed_score)) > 0,
            all(names(final_fit$cell_imputed_score) == names(simulation_res$lineage_assignment)))
  
  metadata <- data.frame(assigned_lineage = simulation_res$lineage_assignment,
                         fatepotential = final_fit$cell_imputed_score)
  rownames(metadata) <- rownames(simulation_res$embedding_mat)
  
  count_mat <- Matrix::Matrix(0, 
                              ncol = nrow(simulation_res$embedding_mat),
                              nrow = 10)
  colnames(count_mat) <- rownames(simulation_res$embedding_mat)
  rownames(count_mat) <- paste0("Gene", 1:10)
  
  seurat_obj <- Seurat::CreateSeuratObject(counts = count_mat, 
                                           meta.data = metadata)
  
  seurat_obj[["fasttopic"]] <- Seurat::CreateDimReducObject(simulation_res$embedding_mat)
  
  seurat_obj@misc$lineage_imputed_count <- final_fit$lineage_imputed_count
  seurat_obj@misc$lineage_observed_count <- simulation_res$lineage_future_size
  
  tmp <- c(simulation_res$coefficient_intercept, simulation_res$coefficient_vec)
  names(tmp) <- c("Intercept", colnames(simulation_res$embedding_mat))
  seurat_obj@misc$true_coefficient_vec <- tmp
  seurat_obj@misc$est_coefficient_vec <- final_fit$coefficient_vec
  
  seurat_obj
}

####################


