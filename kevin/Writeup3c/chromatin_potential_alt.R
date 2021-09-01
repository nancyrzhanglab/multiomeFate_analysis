chromatin_potential_alt <- function(prep_obj, verbose = T,
                                    filepath = NA){
  # pull the appropriate objects for convenience
  mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
  df_x <- prep_obj$df_x; df_y <- prep_obj$df_y
  df_res <- prep_obj$df_res
  snn <- prep_obj$snn
  diffusion_dist <- prep_obj$diffusion_dist
  
  # initialize
  n <- nrow(mat_x)
  
  #: we'll record all the matches and weight-information
  matches_df <- .initialize_matches(mat_x, 
                                    mat_y, 
                                    df_res,
                                    snn,
                                    diffusion_dist)
  iter <- 1
  
  if(!is.na(filepath)) {
    save(mat_x, 
         mat_y, 
         df_res, 
         matches_df, 
         file = filepath)
  }
  
  # while:
  while(any(is.na(df_res$order_rec))){
    if(verbose) print(paste0("Iteration ", iter, ": Recruited percentage (", 
                             round(sum(!is.na(df_res$order_rec))/nrow(df_res), 2), "), Total: ",
                             sum(!is.na(df_res$order_rec)), " cells"))
    if(verbose & !any(is.na(weights))) print(paste0("Weights range from ", round(min(weights),2), " to ", round(max(weights),2)))
    ## estimate res_g
    res <- .estimate_g2(mat_x, 
                        mat_y, 
                        matches_df)
    res_g <- res$res_g
    
    ## construct candidate set
    if(verbose) print("Constructing candidate set")
    res_cand <- .candidate_set2(mat_x, 
                                mat_y, 
                                df_res, 
                                snn,
                                matches_df)
    df_res <- multiomeFate:::.update_chrom_df_cand(df_res, 
                                                   res_cand$vec_cand)
    stopifnot(all(is.na(df_res$order_rec[res_cand$vec_cand])))
    
    ## recruit an element from the candidate se
    if(verbose) print("Recruiting cells")
    enforce_matched <- length(which(df_res$order_rec == 0)) > length(which(df_res$order_rec > 0)) & !bool_oracle
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
    if(verbose) print("Updating matrices")
    matches_df <- .update_matches_df(matches_df, 
                                     res_rec)
    df_res <- .update_chrom_df_rec2(df_res, 
                                    res_rec, 
                                    iter)
    
    if(!is.na(filepath)) {
      save(mat_x, 
           mat_y, 
           df_res, 
           matches_df, 
           res_g,
           file = filepath)
    }
    
    iter <- iter+1
  }
  
  # output
  structure(list(matches_df = matches_df,
                 res_g = res_g, 
                 mat_x = mat_x, 
                 mat_y = mat_y, 
                 df_x = df_x, 
                 df_y = df_y, 
                 df_res = df_res, 
                 snn = snn,
                 diffusion_dist = diffusion_dist,
                 options = options),
            class = "chromatin_potential")
}