.candidate_set2 <- function(mat_x, 
                            mat_y, 
                            df_res, 
                            snn,
                            matches_df){
  # enumerate all metacells that are involved in a match
  matched_idx <- unique(as.numeric(as.matrix(matches_df[which(!matches_df[,"initial_bool"]), 
                                              c("tail", "head")])))
  free_idx <- which(is.na(df_res$order_rec))
  nn_idx <- unique(c(unlist(lapply(matched_idx, function(i){
    which(snn[i,] != 0)
  })), matched_idx))
  
  # find all neighboring metacells that not yet been matched
  list(vec_cand = intersect(nn_idx, free_idx))
}
