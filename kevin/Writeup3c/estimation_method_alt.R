.estimate_g2 <- function(mat_x, 
                         mat_y, 
                         ht_map,
                         matches_df){
  # prepare est_options
  est_options <- .dummy_est_options()
  est_options$ht_map <- ht_map
  
  # form the corresponding regression's X and Y
  tmp <- .form_regression_mat(mat_x, mat_y, matches_df)
  mat_x1 <- tmp$mat_x1
  mat_y2 <- tmp$mat_y2
  
  if("weight" %in% matches_df) weights <- matches_df[,"weight"]
  
  # use the multiomeFate function
  res_g <- multiomeFate:::.estimate_g_glmnet(mat_x1, 
                                             mat_y2, 
                                             weights, 
                                             est_options)
  
  list(res_g = res_g, 
       est_options = est_options)
}

#########

.dummy_est_options <- function(){
  list(method = "threshold_glmnet",
       family = "gaussian", 
       enforce_cis = T, 
       cis_window = 200,
       switch = F, 
       switch_cutoff = 10,
       alpha = 1, 
       standardize = F, 
       intercept = T,
       cv = T, 
       nfolds = 5, 
       initial_quantile = 0.25,
       run_diagnostic = F, 
       hold_initial = F, 
       parallel = F,
       num_iterations = 4,
       cv_choice = "lambda.min",
       verbose = T)
}


.form_regression_mat <- function(mat_x, mat_y, matches_df){
  stopifnot(c("tail", "head") %in% colnames(matches_df))
  
  mat_x1 <- mat_x[matches_df[,"tail"],,drop = F]
  mat_y1 <- mat_y[matches_df[,"tail"],,drop = F]
  mat_y2 <- mat_y[matches_df[,"head"],,drop = F]
  
  list(mat_x1 = mat_x1, 
       mat_y1 = mat_y1, 
       mat_y2 = mat_y2)
}