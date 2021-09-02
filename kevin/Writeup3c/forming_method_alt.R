.initialize_matches <- function(mat_x, 
                                mat_y, 
                                ht_map,
                                df_res,
                                snn,
                                diffusion_dist){
  # find all the initial and terminal states
  r <- max(df_res$init_state, na.rm = T)
  initial_vec <- which(df_res$init_state == -1)
  terminal_list <- lapply(1:r, function(k){
    which(df_res$init_state == k)
  })
  steady_list <- terminal_list
  steady_list[[length(steady_list)+1]] <- initial_vec
  names(steady_list) <- as.character(c(1:r, "-1"))
  steady_vec <- c(initial_vec, unlist(terminal_list))
  
  # for each terminal states, 
  # rank each metacells based on the diffusion distance to all the other 
  # metacells in other initial&termainal states. 
  # Then take the mean among the rankings.
  # Similarly, each initial state, compute the 
  # mean-rank each metacell based on diffusion distance to
  # all the terminal states.
  # Here, smaller rank (closer to 0) means it's closer to all the other cells
  ranking_list <- lapply(1:length(steady_list), function(k){
    idx <- steady_list[[k]]
    
    if(names(steady_list)[k] == "-1"){
      other_vec <- initial_vec
    } else {
      other_vec <- setdiff(steady_vec, idx)
    }
    
    ranking_mat <- diffusion_dist[idx, other_vec, drop = F]
    
    prob_vec <- rep(0, length(other_vec))
    uniq_steady <- unique(df_res[other_vec,"init_state"])
    for(k2 in uniq_steady){
      tmp <- which(df_res[other_vec,"init_state"] == k2)
      prob_vec[tmp] <- 1/length(tmp)
    }
    
    ranking_mat <- diffusion_dist[idx, other_vec, drop = F]
    for(j in 1:ncol(ranking_mat)){
      ranking_mat[,j] <- rank(ranking_mat[,j])
    }
    for(j in 1:nrow(ranking_mat)){
      ranking_mat[j,] <- prob_vec * ranking_mat[j,]
    }
    ranking_vec <- rank(apply(ranking_mat, 1, sum))
    
    data.frame(idx = idx, 
               rank = ranking_vec)
  })
  names(ranking_list) <- names(steady_list)
  
  # Look at snn -- form arrows based on which neighboring
  # metacell has a higher mean-ranking, and if the metacall has
  # the highest rank among its neighbors, it points to itself
  matches_mat <- do.call(rbind, lapply(1:length(steady_list), function(k){
    idx <- ranking_list[[k]][,"idx"]
    do.call(rbind, (lapply(idx, function(i){
      neigh_vec <- intersect(idx, which(snn[i,] != 0))
      neigh_vec <- setdiff(neigh_vec, i)
      
      cell_rank <- ranking_list[[k]][which(idx == i), "rank"]
      neigh_rank <- sapply(neigh_vec, function(j){
        ranking_list[[k]][which(idx == j), "rank"]
      })
      
      if(length(neigh_rank) > 0){
        if(names(ranking_list)[k] == "-1"){
          # if initial, we want cells to point to cells that are closer to other cells
          neigh_vec <- neigh_vec[neigh_rank <= cell_rank]
        } else {
          # if initial, we want cells to point to cells that are further to other cells
          neigh_vec <- neigh_vec[neigh_rank >= cell_rank]
        }
      }
      
      if(length(setdiff(neigh_vec, i)) > 0){
        cbind(i, setdiff(neigh_vec, i))
      } else {
        if(names(ranking_list)[k] == "-1"){
          # do not point initial metacell to itself
          numeric(0)
        } else {
          # point terminal metacell to itself
          c(i, i)
        }
      }
    })))
  }))
  colnames(matches_mat) <- c("tail", "head")
  weight_vec <- sapply(1:nrow(matches_mat), function(i){
    nn <- length(which(snn[matches_mat[i,"tail"],] != 0))
    1/nn
  })
  matches_mat <- cbind(matches_mat, weight_vec)
  colnames(matches_mat)[3] <- "weight"
  # run a check
  stopifnot(length(unique(as.numeric(matches_mat[,c("tail", "head")]))) == length(steady_vec))
  
  # do a first-round of regression
  # (we won't actually return this regression -- its just to get weights)
  tmp <- .form_regression_mat(mat_x, mat_y, matches_mat)
  mat_x1 <- tmp$mat_x1; mat_y1 <- tmp$mat_y1; mat_y2 <- tmp$mat_y2
  res <- .estimate_g2(mat_x1, 
                      mat_y2, 
                      ht_map,
                      matches_mat)
  
  # extract the correlations
  pred_mat <- .predict_yfromx2(mat_x1, res_g)
  cor_vec <- sapply(1:nrow(mat_y2), function(i){
    residual_vec <- mat_y2[i,] - mat_y1[i,]
    if(sum(abs(residual_vec)) <= 1e-6) residual_vec <- stats::runif(length(residual_vec))
    stats::cor(pred_mat[i,] - mat_y1[i,], 
               residual_vec, 
               method = "spearman")
  })
  
  # overwrite the weights -- recompute the weights
  uniq_tail <- unique(matches_mat[,"tail"])
  matches_mat <- do.call(rbind, lapply(uniq_tail, function(i){
    row_idx <- which(matches_mat[,"tail"] == i)
    head_idx <- matches_mat[row_idx, "head"]
    nn <- length(which(snn[i,] != 0))
    tmp_cor <- cor_vec[row_idx]
    
    tmp_cor <- .correlation_to_weight(tmp_cor)
    weight_vec <- tmp_cor/sum(tmp_cor)
    weight_vec <- weight_vec * (1/nn)
    weight_vec
    
    tmp <- cbind(matches_mat[row_idx,c("tail", "head"), drop = F], 
                 cor_vec[row_idx],
                 tmp_cor,
                 weight_vec)
    colnames(tmp) <- c("tail", 
                       "head", 
                       "correlation", 
                       "exp_cor",
                       "weight")
    tmp
  }))
  
  # form the 6-column matrix: tail, head, correlation,
  # exp-cor-weight (using exp(2*x)), regression-weight 
  # (which is exp-cor-weight*1/(nn), used for just the regression),
  # and initial-boolean (which is used only to not count the 
  # initial states as already being matched in df_res, and that
  # the weights learned here will be overwritten in by the "actual"
  # matches in the postprocessing)
  matches_df <- data.frame(tail = matches_mat[,"tail"],
                           head = matches_mat[,"head"],
                           correlation = matches_mat[,"correlation"],
                           exp_cor = matches_mat[,"exp_cor"],
                           weight = matches_mat[,"weight"])
  initial_bool <- sapply(1:nrow(matches_df), function(i){
    any(matches_df[i,c("tail", "head")] %in% initial_vec)
  })
  matches_df[,"initial_bool"] <- initial_bool
  
  matches_df
}

.update_matches_df <- function(matches_df,
                               res_rec){
  # simply append the rows
  rbind(matches_df, res_rec)
}

###########3

.correlation_to_weight <- function(vec){
  exp(2*vec)
}
