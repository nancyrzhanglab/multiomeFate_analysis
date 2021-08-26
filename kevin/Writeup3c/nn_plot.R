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