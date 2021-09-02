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
  list(est_family = "gaussian", 
       est_enforce_cis = T, 
       est_cis_window = 200,
       est_switch = F, 
       est_switch_cutoff = 10,
       est_alpha = 1, 
       est_standardize = F, 
       est_intercept = T,
       est_cv = T, 
       est_nfolds = 5, 
       est_initial_quantile = 0.25,
       est_run_diagnostic = F, 
       est_hold_initial = F, 
       est_parallel = F,
       est_num_iterations = 4,
       est_cv_choice = "lambda.min",
       est_verbose = T)
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