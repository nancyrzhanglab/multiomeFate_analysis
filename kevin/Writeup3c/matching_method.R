.make_matches <- function(mat_x, 
                          mat_y, 
                          vec_cand, 
                          res_g,
                          df_res,
                          snn,
                          diffusion_dist,
                          gene_weights,
                          factor = 1.5){
  len <- length(vec_cand)
  
  # for each candidate metacell, look at its diffusion distance to 
  # terminal/initial states
  uniq_steady <- unique(df_res[!is.na(df_res[,"init_state"]), "init_state"])
  steady_list <- lapply(uniq_steady, function(i){
    which(df_res[,"init_state"] == i)
  })
  names(steady_list) <- uniq_steady
  dist_to_steady <- sapply(steady_list, function(steady_vec){
    matrixStats::rowMeans2(diffusion_dist[vec_cand, steady_vec, drop = F])
  })
  rownames(dist_to_steady) <- vec_cand
  colnames(dist_to_steady) <- uniq_steady
  close_mat <- t(sapply(1:len, function(i){
    idx <- which.min(dist_to_steady[i,])
    bool <- dist_to_steady[i,idx] <= min(dist_to_steady[i,-idx])/factor
    if(bool){
      c(1, uniq_steady[idx])
    } else {
      c(0, NA)
    }
  }))
  close_df <- data.frame(bool = as.logical(close_mat[,1]),
                         steady_label = close_mat[,2])
  
  # compute this metacell's predicted RNA expression
  # and compute the change in RNA expression
  pred_mat <- .predict_yfromx2(mat_x[vec_cand,,drop = F], 
                               res_g)
  
  potential_list <- lapply(1:len, function(i){
    idx <- vec_cand[i]
    nn <- which(snn[idx,] != 0)
    
    # if it's closer to one terminal state than any other terminal
    # state or initial state (say, by a factor of 1.5),
    # then only consider neighboring metacells with
    # closer diffusion distance to said terminal state
    if(close_df[i,"bool"] && (is.na(df_res[idx,"init_state"]) || df_res[idx,"init_state"] != "-1")){
      steady_label <- close_df[i,"steady_label"]
      steady_vec <- steady_list[[which(names(steady_list) == steady_label)]]
      steady_dist <- matrixStats::rowMeans2(diffusion_dist[nn, steady_vec, drop = F])
      
      if(steady_label != "-1"){
        if(any(steady_dist < min(dist_to_steady[i,]))){
          nn <- nn[which(steady_dist < min(dist_to_steady[i,]))]
        } 
      } else {
        # OTHERWISE 1: if the metacell is closer to the initial 
        # state than any other terminal state (using the same heuristic),
        # it can only point to cells further way from the initial
        # state in diffusion-distance
        if(any(steady_dist > min(dist_to_steady[i,]))){
          nn <- nn[steady_dist > min(dist_to_steady[i,])]
        } 
      }
      
      # OTHERWISE 2: if the metacell is an initial state, can
      # only point to neighboring metacells that are 1) closer to
      # the terminal states based on diffusion distance or 2)
      # cells that are not an initial state
    } else if(!is.na(df_res[idx,"init_state"]) && df_res[idx,"init_state"] == "-1"){
      dist_vec <- sapply(steady_list, function(steady_vec){
        matrixStats::rowMeans2(diffusion_dist[nn, steady_vec, drop = F])
      })
      
      if(any(dist_vec < min(dist_to_steady[i,]))){
        nn <- nn[dist_vec < min(dist_to_steady[i,])]
      }
      
      # OTHERWISE 3:, consider all neighboring metacells that 
      # are NOT in the initial state
    } else {
      init_idx <- which(df_res[nn,"init_state"] == "-1")
      if(length(init_idx) > 0 & length(init_idx) < length(nn)){
        nn <- nn[-init_idx]
      }
    }
    
    as.numeric(nn)
  })
  stopifnot(all(sapply(potential_list, length) > 0))
  
  # among all considered neighboring metacells, compare their
  # RNA expression - this metacell's RNA
  cor_list <- lapply(1:len, function(i){
    residual_vec1 <- pred_mat[i,] - mat_y[vec_cand[i],,drop = F]
    nn_idx <- potential_list[[i]]
    sapply(nn_idx, function(j){
      residual_vec2 <- mat_y[j,] - mat_y[vec_cand[i],,drop = F]
      
      # compute the correlation among the difference vectors
      wCorr::weightedCorr(as.numeric(residual_vec1), 
                         as.numeric(residual_vec2),
                         weights = gene_weights,
                         method = "Pearson")
    })
  })
  
  # form the appropriate weights
  weights <- vector("list", length = len)
  for(i in 1:len){
    idx <- which(cor_list[[i]] < 0)
    if(length(idx) > 0){
      potential_list[[i]] <- potential_list[[i]][-idx]
      cor_list[[i]] <- cor_list[[i]][-idx]
    }
    
    nn <- length(which(snn[vec_cand[i],] != 0))
    tmp <- .correlation_to_weight(cor_list[[i]])
    weights[[i]] <- (tmp/sum(tmp))*(1/nn)
  }
  
  # return the appropriate additional rows to be appended to matched_mat
  matches_mat <- do.call(rbind, lapply(1:len, function(i){
    cbind(rep(vec_cand[i], length(potential_list[[i]])),
          potential_list[[i]],
          cor_list[[i]],
          .correlation_to_weight(cor_list[[i]]),
          weights[[i]])
  }))
  matches_df <- as.data.frame(matches_mat)
  colnames(matches_df) <- c("tail",
                            "head",
                            "correlation",
                            "exp_cor",
                            "weight")
  matches_df[,"initial_bool"] <- rep(F, nrow(matches_df))
  
  matches_df
}


.update_chrom_df_rec2 <- function(df_res, 
                                  res_rec, 
                                  iter){
  # mark all metacells that got matched
  idx <- which(is.na(df_res[,"order_rec"]))
  idx <- intersect(idx, res_rec[,"tail"])
  df_res[idx,"order_rec"] <- iter
  
  df_res
}
