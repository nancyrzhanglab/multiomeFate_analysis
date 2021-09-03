chromatin_potential_postprocess <- function(chrom_obj){
  matches_df <- chrom_obj$matches_df
  df_res <- chrom_obj$df_res
  snn <- chrom_obj$snn
  diffusion_dist <- chrom_obj$diffusion_dist
  
  # form the metacell weighted graph where the weights are
  # the exp-correlation-proportions.
  # for any missing edges that do not involve any initial/terminal
  # cells, put exp(-2) [the minimum possible weight]
  n <- nrow(snn)
  adj_mat <- matrix(0, n, n)
  adj_mat[adj_mat != 0] <- exp(-2)
  for(row_idx in 1:nrow(matches_df)){
    i <- matches_df[row_idx,"tail"]
    j <- matches_df[row_idx,"head"]
    adj_mat[i,j] <- matches_df[row_idx,"exp_cor"]
  }
  for(i in 1:n){
    if(!is.na(df_res[i,"init_state"]) && df_res[i,"init_state"] != "-1") next()
    
    idx1 <- which(adj_mat[i,] == 0)
    idx2 <- which(snn[i,] != 0)
    adj_mat[i, intersect(idx1, idx2)] <- exp(-2)
  }
  
  for(i in 1:n){
    adj_mat[i,] <- adj_mat[i,]/sum(adj_mat[i,])
  }
  
  r <- max(df_res$init_state, na.rm = T)
  terminal_list <- lapply(1:r, function(k){
    which(df_res$init_state == k)
  })
  terminal_vec <- unlist(terminal_list)
  
  A <- solve(diag(n - length(terminal_vec)) - adj_mat[-terminal_vec, -terminal_vec]) %*% adj_mat[-terminal_vec, terminal_vec]
  colnames(A) <- terminal_vec
  rownames(A) <- colnames(diffusion_dist)[-terminal_vec]
  
  fate_prob <- sapply(1:r, function(k){
    idx <- which(colnames(A) %in% as.character(terminal_list[[k]]))
    rowSums(A[,idx,drop = F])
  })
  colnames(fate_prob) <- as.character(1:r)
  
  # [[can we manually set the initial state to be near-uniform?]]
  # [[... this sounds super janky]]
  
  fate_prob
}