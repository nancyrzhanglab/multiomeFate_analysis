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