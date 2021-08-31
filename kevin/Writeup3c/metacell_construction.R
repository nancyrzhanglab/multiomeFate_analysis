form_metacell_matrix <- function(dat, clustering){
  stopifnot(is.character(clustering), length(clustering) == nrow(dat))
  
  uniq_clust <- sort(unique(clustering))
  clust_mat <- t(sapply(uniq_clust, function(clust){
    idx <- which(clustering == clust)
    apply(dat[idx,,drop = F], 2, median)
  }))
  rownames(clust_mat) <- uniq_clust
  
  clust_mat
}

form_snn_graph <- function(dat, k, distance){
  n <- nrow(dat)
  
  set.seed(10)
  nn_res <- knn.covertree::find_knn(dat, 
                                    k = 10,
                                    distance = "cosine")

  adj_mat <- matrix(0, n, n)
  for(i in 1:n){
    adj_mat[i, nn_res$index[i,]] <- 1
  }
  diag(adj_mat) <- 0
  adj_mat <- adj_mat*t(adj_mat)
  
  idx <- which(rowSums(adj_mat) == 0)
  for(i in idx){
    adj_mat[i,nn_res$nn.idx[i,2]] <- 1
  }
  
  adj_mat <- adj_mat + t(adj_mat)
  adj_mat[adj_mat > 0] <- 1
  diag(adj_mat) <- 0
  
  rownames(adj_mat) <- rownames(dat)
  colnames(adj_mat) <- rownames(dat)
  
  adj_mat
}
