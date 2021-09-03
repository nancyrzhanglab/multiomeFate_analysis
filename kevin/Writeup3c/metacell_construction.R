form_metacell_matrix <- function(dat, clustering, func = median){
  stopifnot(is.character(clustering), length(clustering) == nrow(dat))
  
  uniq_clust <- sort(unique(clustering))
  clust_mat <- t(sapply(uniq_clust, function(clust){
    idx <- which(clustering == clust)
    apply(dat[idx,,drop = F], 2, func)
  }))
  rownames(clust_mat) <- uniq_clust
  
  clust_mat
}

form_snn_graph <- function(dat, initial_vec, terminal_list){
  stopifnot(is.character(initial_vec), all(sapply(terminal_list, is.character)))
  n <- nrow(dat)
  
  dist_mat <- acos(lsa::cosine(t(dat)))/pi
  rownames(dist_mat) <- NULL
  colnames(dist_mat) <- NULL
  
  # create weighted graph
  g <- igraph::graph_from_adjacency_matrix(dist_mat,
                                           mode = "undirected",
                                           weighted = T)
  g2 <- igraph::mst(g, algorithm = "prim")
  # ensures connectivity
  mst_edge_mat <- igraph::as_edgelist(g2)
  
  # ensures enough edges to start with
  nn_res2 <- knn.covertree::find_knn(dat, 
                                     k = 3,
                                     distance = "cosine")
  
  # now ensure that the terminal states are connected
  k <- 3
  while(TRUE){
    set.seed(10)
    nn_res <- knn.covertree::find_knn(dat, 
                                      k = k,
                                      distance = "cosine")
    
    adj_mat <- .form_snn_from_edgelists(n, 
                                        mst_edge_mat, 
                                        nn_res2$index,
                                        nn_res$index,
                                        rownames(dat))
    
    if(.check_snn(adj_mat, initial_vec, terminal_list)) break()
    
    k <- k+1
  }
 
  snn <- adj_mat
  for(i in 1:n){
    idx <- which(adj_mat[i,] != 0)
    snn[i,idx] <- dist_mat[i,idx]
  }
  
  rownames(snn) <- rownames(dat)
  colnames(snn) <- rownames(dat)
  
  # snn has non-zero values that represent the distance
  list(snn = snn, adj_mat = adj_mat)
}

###########################

.form_snn_from_edgelists <- function(n, 
                                     mst_edge_mat, 
                                     knn2_edge_mat,
                                     knn_edge_mat,
                                     rowname_vec){
  adj_mat1 <- matrix(0, n, n)
  for(i in 1:nrow(mst_edge_mat)){
    idx1 <- mst_edge_mat[i,1]
    idx2 <- mst_edge_mat[i,2]
    
    adj_mat1[idx1, idx2] <- 1
    adj_mat1[idx2, idx1] <- 1
  }
  
  adj_mat2 <- matrix(0, n, n)
  for(i in 1:nrow(knn2_edge_mat)){
    adj_mat2[i, knn2_edge_mat[i,]] <- 1
  }
  diag(adj_mat2) <- 0
  adj_mat2 <- adj_mat2 + t(adj_mat2)
  
  adj_mat3 <- matrix(0, n, n)
  for(i in 1:nrow(knn_edge_mat)){
    adj_mat3[i, knn_edge_mat[i,]] <- 1
  }
  diag(adj_mat3) <- 0
  adj_mat3 <- adj_mat3*t(adj_mat3)
  
  adj_mat <- adj_mat1 + adj_mat2 + adj_mat3
  adj_mat[adj_mat != 0] <- 1
  
  rownames(adj_mat) <- rowname_vec
  colnames(adj_mat) <- rowname_vec
  
  adj_mat
}


.check_snn <- function(adj_mat, 
                       initial_vec, 
                       terminal_list,
                       check_steady = T){
  steady_vec <- c(initial_vec, unlist(terminal_list))
  
  # check that only steady states have degree 1
  col_vec <- colSums(adj_mat)
  bool1 <- all(col_vec[-which(colnames(adj_mat) %in% steady_vec)] > 1)
  
  if(check_steady){
    steady_list <- c(list(initial_vec), terminal_list)
    adj_mat2 <- adj_mat
    rownames(adj_mat2) <- NULL
    colnames(adj_mat2) <- NULL
    g <- igraph::graph_from_adjacency_matrix(adj_mat2, 
                                             mode = "undirected")
    bool_vec <- sapply(1:length(steady_list), function(i){
      if(length(steady_list[[i]]) == 1) return(TRUE)
      g2 <- igraph::induced_subgraph(g, which(rownames(adj_mat) %in% steady_list[[i]]))
      igraph::components(g2)$no == 1
    })
    
    return(all(c(bool1, bool_vec)))
  } else {
    return(bool1)
  }
}
