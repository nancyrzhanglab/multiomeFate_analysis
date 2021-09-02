mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
df_x <- prep_obj$df_x; df_y <- prep_obj$df_y
ht_map <- prep_obj$ht_map
df_res <- prep_obj$df_res
snn <- prep_obj$snn
diffusion_dist <- prep_obj$diffusion_dist

# initialize
n <- nrow(mat_x)

#: we'll record all the matches and weight-information
matches_df <- .initialize_matches(mat_x,
                                  mat_y,
                                  ht_map,
                                  df_res,
                                  snn,
                                  diffusion_dist)

iter <- 1

# while:
while(any(is.na(df_res$order_rec))){
  print(paste0("Iteration ", iter, ": Recruited percentage (", 
               round(sum(!is.na(df_res$order_rec))/nrow(df_res), 2), "), Total: ",
               sum(!is.na(df_res$order_rec)), " cells"))
  ## estimate res_g
  res_g <- .estimate_g2(mat_x, 
                        mat_y, 
                        ht_map,
                        matches_df)
  
  ## construct candidate set
  print("Constructing candidate set")
  res_cand <- .candidate_set2(mat_x, 
                              mat_y, 
                              df_res, 
                              snn,
                              matches_df)
  df_res <- multiomeFate:::.update_chrom_df_cand(df_res, 
                                                 res_cand$vec_cand)
  stopifnot(all(is.na(df_res$order_rec[res_cand$vec_cand])))
  
  ## recruit an element from the candidate se
  print("Recruiting cells")
  #[NOTE: Rebrand this step as matching]
  #[NOTE: Is it possible that cell-level information would help with this step?]
  res_rec <- .make_matches(mat_x, 
                           mat_y, 
                           res_cand$vec_cand,
                           res_g,
                           df_res,
                           snn,
                           diffusion_dist)
  
  ## update
  print("Updating matrices")
  matches_df <- .update_matches_df(matches_df, 
                                   res_rec)
  df_res <- .update_chrom_df_rec2(df_res, 
                                  res_rec, 
                                  iter)
  
  iter <- iter+1
}