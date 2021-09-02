# [[DISCONTINUED -- since I can't figure out how to get
# this to work...]]

# I don't use directed diffusion distance since it doesn't
# guarantee that the initial states are the most-uniform

# we might want to use Markov-logic here eventually?
# say: smooth the estimated markov-transition matrix via rank-K
# and then use absorption probabilities? 
chromatin_potential_postprocess <- function(chrom_obj, max_iter = 100){
  matches_df <- chrom_obj$matches_df
  df_res <- chrom_obj$df_res
  snn <- chrom_obj$snn
  diffusion_dist <- chrom_obj$diffusion_dist
  
  # form the metacell weighted graph where the weights are
  # the exp-correlation-proportions.
  # for any missing edges that do not involve any initial/terminal
  # cells, put exp(-2) [the minimum possible weight]
  n <- nrow(snn)
  adj_mat <- snn
  adj_mat[adj_mat != 0] <- exp(-2)
  for(row_idx in 1:nrow(matches_df)){
    i <- matches_df[row_idx,"tail"]
    j <- matches_df[row_idx,"head"]
    adj_mat[i,j] <- matches_df[row_idx,"exp_cor"]
  }
  
  # under the assumption that cells should be more differented
  # as time proceeds, use the following logic:
  # set the fate-probability vector of all terminal states
  # to be an indicator vector, and the initial state to be
  # 1/K vector, and all other cells NA
  r <- max(df_res$init_state, na.rm = T)
  fate_prob <- matrix(NA, n, r)
  initial_vec <- which(df_res$init_state == -1)
  terminal_list <- lapply(1:r, function(k){
    which(df_res$init_state == k)
  })
  steady_vec <- c(initial_vec, unlist(terminal_list))
  for(i in initial_vec){
    fate_prob[i,] <- rep(1/r, r)
  }
  for(k in 1:length(terminal_list)){
    for(i in terminal_list[[k]]){
      fate_prob[i,-k] <- 0
      fate_prob[i,k] <- 1
    }
  }
  
  # order all the cells based on the min-distance to any initial/terminal state
  # this is the order we'll deal with cells
  cell_order <- order(apply(diffusion_dist, 1, function(x){
    min(x[steady_vec])
  }), decreasing = F)
  cell_order <- setdiff(cell_order, steady_vec)

  
  # this proceeds using in the following steps:
  
  # using the posited ordering, initialize via backpropogating, where we 
  # slowly flip the NAs
  # and set each metacell equal to the weighted sum of its outarrows
  # (this is an iterative process)
  # we go in /reverse/ here to avoid vectors from being repeated all the time
  while(any(is.na(fate_prob))){
    for(i in rev(cell_order)){
      if(!is.na(fate_prob[i,1])) next()
      nn <- which(adj_mat[i,] != 0)
      if(all(is.na(fate_prob[nn,1]))) next()
      
      val <- adj_mat[i,nn[!is.na(fate_prob[nn,1])]]
      if(length(val) == 1) diag_mat <- matrix(val, 1, 1) else diag_mat <- diag(val)
      nn <- nn[!is.na(fate_prob[nn,1])]
      fate_prob[i,] <- colSums(diag_mat %*% fate_prob[nn,,drop = F])/sum(val)
    }
  }
  
  iter <- 1
  while(TRUE){
    bool_vec <- rep(FALSE, length(cell_order))
    
    # then keep iterating under the two rules:
    # we are trying to find a vector for each metacell such that:
    # 1) each metacell's vector is more-uniform than the
    # weighted sum of its out-arrows
    # 2) each metacell is less-uniform than the weighted sum
    # of its in-arrows
    for(i in 1:length(cell_order)){
      cell_idx <- cell_order[i]
      
      nn <- which(adj_mat[cell_idx,] != 0)
      val <- adj_mat[cell_idx,nn]
      if(length(val) == 1) diag_mat <- matrix(val, 1, 1) else diag_mat <- diag(val)
      out_vec <- colSums(diag_mat %*% fate_prob[nn,,drop = F])/sum(val)
      
      nn <- which(adj_mat[,cell_idx] != 0)
      val <- adj_mat[nn,cell_idx]
      if(length(val) == 1) diag_mat <- matrix(val, 1, 1) else diag_mat <- diag(val)
      in_vec <- colSums(diag_mat %*% fate_prob[nn,,drop = F])/sum(val)
      
      tmp <- apply(rbind(out_vec, in_vec), 2, range)
      bool1 <- all(sapply(1:r, function(j){
        fate_prob[cell_idx,j] >= tmp[1,j] & fate_prob[cell_idx,j] <= tmp[2,j]
      }))
      cell_uniformity <- .uniformity(fate_prob[cell_idx,])
      in_uniformity <- .uniformity(in_vec)
      out_uniformity <- .uniformity(out_vec)
      bool2 <- cell_uniformity <= out_uniformity & cell_uniformity >= out_uniformity
      
      bool_vec[i] <- bool1 & bool2
      
      if(!bool_vec[i]){
        # then, based on the iteration, repeat the following steps:
        # MODIFYING
        # using the posited ordering, for each compute the weighted
        # average of its in-arrows and average of its out-arrows
        # this gives 2 K-dimension vectors, and the cell has its own K-dimensional
        # vector
        # perform some optimization so that within the values
        # set by the 2-K dimension vectors AND is sandwiched in terms of
        # it's non-uniformity
        
      }
    }
    
    print(paste0("On iteration ", iter, " with ", length(which(!bool_vec)), " violations"))
    
    # TERMINATION?
    # see if our fate-prob vectors satsifies our system-of-equations
    # if not, repeat the above process again
    
    if(all(bool_vec) | iter > max_iter) break()
    iter <- iter + 1
  }
}

.optimization_uniformity <- function(vec, 
                                     lower_range, 
                                     upper_range,
                                     in_uniformity,
                                     out_uniformity,
                                     tol = 1e-6){
  stopifnot(out_uniformity + tol >= in_uniformity)
  
  r <- length(vec)
  obj_func <- optiSolve::quadfun(Q = diag(r), 
                                 a = -2*vec)
  lower_con <- optiSolve::lbcon(lower_range)
  upper_con <- optiSolve::ubcon(upper_range)
  
  Q <- 
  out_quad_con <- optiSolve::quadcon(Q = diag(r),
                                     a = rep(-2/r, r),
                                     d = 1/r,
                                     dir = "<=",
                                     val = out_uniformity)
  in_quad_con <- optiSolve::quadcon(Q = -diag(r),
                                    a = rep(2/r, r),
                                    d = -1/r,
                                    dir = "<=",
                                    val = -in_uniformity)
}


# smaller number equates to more uniform
.uniformity <- function(vec){
  r <- length(vec)
  sum((vec - 1/r)^2)
}