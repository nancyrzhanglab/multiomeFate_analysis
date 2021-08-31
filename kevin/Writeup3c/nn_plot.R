nn_plot <- function(seurat_obj, nn_mat, celltype, end_states){
  stopifnot(length(end_states) == 4, length(names(end_states)) == 4)
  
  n <- nrow(nn_mat)
  mat <- matrix(0, n, n)
  for(i in 1:n){
    mat[i,nn_mat[i,]] <- 1
  }
  
  g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected")
  g <- igraph::simplify(g)
  
  dist_mat <- sapply(end_states, function(x){
    idx <- which(celltype %in% x)
    tmp <- igraph::distances(g, to = idx)
    apply(tmp, 1, mean)
  })
  max_dist <- max(dist_mat)
  
  rownames(dist_mat) <- rownames(seurat_obj@meta.data)
  colnames(dist_mat) <- paste0("nndist_", 1:ncol(dist_mat))
  seurat_obj[["nndist"]] <- Seurat::CreateDimReducObject(embedding = dist_mat, key = "nndist_")
  
  plot_list <- lapply(1:ncol(dist_mat), function(i){
    plot1 <- Seurat::FeaturePlot(seurat_obj, features = paste0("nndist_", i), reduction = "wnn.umap")
    plot1 <- plot1 + ggplot2::labs(title = paste0("NN distance to ", names(end_states)[i])) + 
      ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = c(0,max_dist))
    plot1
  })
  
  tmp2 <- cowplot::plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]])
  
  list(g = g, cowplot = tmp2)
}

dirnn_plot <- function(seurat_obj, list_diagnos, celltype, end_states){
  stopifnot(length(end_states) == 4, length(names(end_states)) == 4)
  n <- nrow(seurat_obj@meta.data)
  
  list_diagnos <- res$list_diagnos
  directed_mat <- matrix(0, n, n)
  for(iter in 1:length(list_diagnos)){
    for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
      idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
      idx_to <- list_diagnos[[iter]]$recruit$postprocess[[i]]$to
      
      directed_mat[idx_from,idx_to] <- 1
    }
  }
  # also account for the terminal cells
  for(end_state in end_states){
    idx <- which(celltype %in% end_state)
    directed_mat[idx,idx] <- 1
  }
  
  g <- igraph::graph_from_adjacency_matrix(directed_mat, mode = "directed")
  
  dist_mat <- sapply(end_states, function(x){
    idx <- which(celltype %in% x)
    tmp <- igraph::distances(g, to = idx, mode = "out", weights = NA)
    apply(tmp, 1, function(x){
      if(all(is.infinite(x))) return(NA) 
      mean(x[!is.infinite(x)])
    })
  })
  max_dist <- max(dist_mat, na.rm = T)
  
  rownames(dist_mat) <- rownames(seurat_obj@meta.data)
  colnames(dist_mat) <- paste0("dirnndist_", 1:ncol(dist_mat))
  seurat_obj[["dirnndist"]] <- Seurat::CreateDimReducObject(embedding = dist_mat, key = "dirnndist_")
  
  plot_list <- lapply(1:ncol(dist_mat), function(i){
    plot1 <- Seurat::FeaturePlot(seurat_obj, features = paste0("dirnndist_", i), reduction = "wnn.umap")
    plot1 <- plot1 + ggplot2::labs(title = paste0("Direct NN dist. to ", names(end_states)[i])) + 
      ggplot2::scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 9, name = "Purples"), limits = c(0,max_dist))
    plot1
  })
  
  tmp2 <- cowplot::plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]])
  
  list(g = g, cowplot = tmp2)
}

compute_nndist <- function(mat_x, mat_y, dim_method = "pca", nn_method = "annoy",
                           form_method = "average", est_method = "glmnet",
                           cand_method = "nn_any", rec_method = "distant_cor",
                           options){
  full_options <- multiomeFate:::.chrom_options(dim_method, nn_method,
                                 form_method, est_method, 
                                 cand_method, rec_method, 
                                 options)
  dim_options <- full_options$dim_options; nn_options <- full_options$nn_options
  
  set.seed(10)
  dim_reduc_obj <- vector("list", 0)
  tmp <- multiomeFate:::dimension_reduction(mat_x, mode = "x", dim_options)
  x_scores <- tmp$scores; dim_reduc_obj$x <- tmp$dim_reduc_obj
  tmp <- multiomeFate:::dimension_reduction(mat_y, mode = "y", dim_options)
  y_scores <- tmp$scores; dim_reduc_obj$y <- tmp$dim_reduc_obj
  
  # form the nn
  n <- nrow(mat_x)
  if(nn_options$include_x & nn_options$include_y){
    all_scores <- cbind(x_scores, y_scores)
  } else if(nn_options$include_x){
    all_scores <- x_scores
  } else{
    all_scores <- y_scores
  }
  
  set.seed(10)
  nn_obj <- multiomeFate:::.nearest_neighbor_annoy(all_scores, nn_options)
  
  vec <- sapply(1:n, function(i){
    tmp <- nn_obj$getNNsByItemList(i-1, nn_options$nn+1, earch_k = -1, include_distances = T)
    max(tmp$distance)
  })
  
  list(all_scores = all_scores, vec = vec)
}

density_plot <- function(seurat_obj, mat_x, mat_y,
                         nn_vec, include_x = T, include_y = T, 
                         nn_metric = "cosine",
                         rank_x = 50, rank_y = 30){
  stopifnot(length(nn_vec) == 4)
  
  dens_mat <- sapply(nn_vec, function(nn){
    compute_nndist(mat_x, mat_y, options = list(nn_nn = nn, nn_metric = nn_metric,
                                                dim_dims_x = 2:rank_x,
                                                dim_dims_y = 1:rank_y, 
                                                nn_include_x = include_x,
                                                nn_include_y = include_y))$vec
  })
  
  
  rownames(dens_mat) <- rownames(seurat_obj@meta.data)
  colnames(dens_mat) <- paste0("dens_", 1:ncol(dens_mat))
  seurat_obj[["dens"]] <- Seurat::CreateDimReducObject(embedding = dens_mat, key = "dens_")
  
  plot_list <- lapply(1:ncol(dens_mat), function(i){
    plot1 <- Seurat::FeaturePlot(seurat_obj, features = paste0("dens_", i), reduction = "wnn.umap")
    plot1 <- plot1 + ggplot2::labs(title = paste0("Density to ", nn_vec[i], " NN"))
    plot1
  })
  
  tmp2 <- cowplot::plot_grid(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]])
  
  list(dens_mat = dens_mat, cowplot = tmp2)
}