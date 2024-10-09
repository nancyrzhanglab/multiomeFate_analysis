.form_simulation_seurat_fate <- function(final_fit,
                                         simulation_res){
  metadata <- cbind(simulation_res$lineage_assignment,
                    final_fit$cell_imputed_score)
  colnames(metadata) <- c("assigned_lineage", "fatepotential")
  rownames(metadata) <- rownames(simulation_res$embedding_mat)
  metadata <- data.frame(metadata)
  metadata$assigned_lineage <- factor(paste0("lineage:", metadata$assigned_lineage)) 
  
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

.plot_anova_helper <- function(seurat_object,
                               cell_imputed_score,
                               assigned_lineage_variable,
                               lineage_future_size,
                               bool_mark_mean = TRUE,
                               bool_mark_median = TRUE,
                               min_lineage_size = 2,
                               num_lineages = 20,
                               ylab = "",
                               ylim = NA){
  stopifnot(length(names(cell_imputed_score)) == length(cell_imputed_score))
  
  if(any(is.na(cell_imputed_score))){
    cell_imputed_score <- cell_imputed_score[!is.na(cell_imputed_score)]
  }
  
  # grab the vector of which celltype-time each cell is
  assigned_lineage <- seurat_object@meta.data[,assigned_lineage_variable]
  names(assigned_lineage) <- Seurat::Cells(seurat_object)
  
  # determine which lineages qualify to be in the plot
  lineage_vec <- assigned_lineage[names(cell_imputed_score)]
  tab_vec <- table(assigned_lineage)
  tab_vec <- tab_vec[tab_vec >= min_lineage_size]
  lineage_names <- names(lineage_future_size)[order(lineage_future_size, decreasing = TRUE)[1:num_lineages]]
  idx <- which(lineage_vec %in% lineage_names)
  
  # form data frame
  df <- data.frame(lineage = lineage_vec[idx],
                   imputed_count = cell_imputed_score[idx])
  df_tmp <- df; df_tmp$lineage <- as.factor(df_tmp$lineage)
  anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    imputed_count = cell_imputed_score)
  df <- rbind(df, df2)
  
  # compute percentage
  lineage_effect <- multiomeFate:::.anova_percentage(
    df = df_tmp,
    lineage_variable = "lineage",
    value_variable = "imputed_count"
  )
  
  col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
  names(col_vec) <- c(lineage_names, "All")
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=imputed_count))
  plot1 <- plot1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=lineage))
  plot1 <- plot1 + ggplot2::scale_fill_manual(values = col_vec) 
  plot1 <- plot1 + ggplot2::geom_jitter(shape=16, 
                                        position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  plot1 <- plot1 + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::geom_boxplot(width=0.05)
  plot1 <- plot1 + ggplot2::scale_x_discrete(limits = c(lineage_names, "All"),
                                             guide = ggplot2::guide_axis(angle = 45))
  plot1 <- plot1 + ggplot2::ylab(ylab)
  
  if(!all(is.na(ylim))){
    plot1 <- plot1 + ggplot2::ylim(ylim[1], ylim[2])
  }
  
  if(bool_mark_mean) 
    plot1 <- plot1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  
  if(bool_mark_median) 
    plot1 <- plot1 + ggplot2::stat_summary(fun=max, geom="point", shape=10, size=5, color="blue")
  
  plot1 <- plot1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2), 
                                           ", Lineage effect = ", lineage_effect, "%"))
  
  plot1
}
