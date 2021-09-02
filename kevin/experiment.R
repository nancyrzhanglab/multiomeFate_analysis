mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
df_x <- prep_obj$df_x; df_y <- prep_obj$df_y
ht_map <- prep_obj$ht_map
df_res <- prep_obj$df_res
snn <- prep_obj$snn
diffusion_dist <- prep_obj$diffusion_dist

# initialize
n <- nrow(mat_x)

#: we'll record all the matches and weight-information
# matches_df <- .initialize_matches(mat_x, 
#                                   mat_y, 
#                                   ht_map,
#                                   df_res,
#                                   snn,
#                                   diffusion_dist)

##############
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

tmp <- .form_regression_mat(mat_x, mat_y, matches_mat)
mat_x1 <- tmp$mat_x1; mat_y1 <- tmp$mat_y1; mat_y2 <- tmp$mat_y2
res <- .estimate_g2(mat_x1, 
                    mat_y2, 
                    ht_map,
                    matches_mat)