mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
df_x <- prep_obj$df_x; df_y <- prep_obj$df_y
ht_map <- prep_obj$ht_map
df_res <- prep_obj$df_res
snn <- prep_obj$snn
diffusion_dist <- prep_obj$diffusion_dist

# initialize
n <- nrow(mat_x)

# #: we'll record all the matches and weight-information
# matches_df <- .initialize_matches(mat_x, 
#                                   mat_y, 
#                                   ht_map,
#                                   df_res,
#                                   snn,
#                                   diffusion_dist)
# iter <- 1

########################

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
  other_vec <- setdiff(steady_vec, idx)
  
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

lapply(ranking_list, function(x){
  x[,1] <- as.numeric(rownames(diffusion_dist)[x[,1]])
  x
})

