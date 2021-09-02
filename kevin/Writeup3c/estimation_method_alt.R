.estimate_g2 <- function(mat_x, 
                         mat_y, 
                         ht_map,
                         matches_df){

  # form the corresponding regression's X and Y
  tmp <- .form_regression_mat(mat_x, mat_y, matches_df)
  mat_x1 <- tmp$mat_x1
  mat_y2 <- tmp$mat_y2
  
  if("weight" %in% matches_df) weights <- matches_df[,"weight"]
  
  # use the multiomeFate function
  res_g <- .estimate_g_glmnet2(mat_x1, 
                               mat_y2, 
                               weights, 
                               ht_map)
  
  list(res_g = res_g, 
       est_options = est_options)
}

#########

.form_regression_mat <- function(mat_x, mat_y, matches_df){
  stopifnot(c("tail", "head") %in% colnames(matches_df))
  
  mat_x1 <- mat_x[matches_df[,"tail"],,drop = F]
  mat_y1 <- mat_y[matches_df[,"tail"],,drop = F]
  mat_y2 <- mat_y[matches_df[,"head"],,drop = F]
  
  list(mat_x1 = mat_x1, 
       mat_y1 = mat_y1, 
       mat_y2 = mat_y2)
}


.estimate_g_glmnet2 <- function(mat_x1, 
                                mat_y2, 
                                weights, 
                                ht_map,
                                verbose = T){
  stopifnot(nrow(mat_x1) == nrow(mat_y2))
  p1 <- ncol(mat_x1); p2 <- ncol(mat_y2)
  
  # initialize variables for the loop
  my_lapply <- pbapply::pblapply
  if(verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
  
  list_res <- my_lapply(1:p2, function(j){
    idx_x <- est_options$ht_map[[as.character(j)]]
    if(length(idx_x) == 0) return(val_int = mean(mat_y2[,j]), 
                                  vec_coef = rep(0, p1))
    idx_y <- j
    
    ## apply glmnet
    .threshold_glmnet_fancy2(mat_x1[,idx_x,drop = F], 
                             mat_y2[,idx_y],
                             weights = weights)
  })
  
  .transform_est_matrix2(list_res, ht_map, p1)
}

##############

.glmnet_fancy2 <- function(x, y, 
                           weights, 
                           tol = 1e-2){
  n <- length(y); p <- ncol(x)
  
  if(nrow(x) == 1 || all(matrixStats::colSds(x) <= tol) || stats::sd(y) <= tol){
    return(list(val_int = mean(y), vec_coef = rep(0, ncol(x))))
  }
  
  # use LM
  df <- as.data.frame(cbind(y, x))
  colnames(df) <- c("y", paste0("x", 1:ncol(x)))
  lm_formula <- stats::as.formula("y ~ .")
  lm_fit <- stats::lm(lm_formula, data = df, weights = weights)
  
  list(val_int = as.numeric(lm$coefficients[1]), 
       vec_coef = as.numeric(lm$coefficients[-1]))
}

.threshold_glmnet_fancy2 <- function(x, y, 
                                    weights, 
                                    tol = 1e-2){
  if(nrow(x) == 1 || all(matrixStats::colSds(x) <= tol) || stats::sd(y) <= tol){
    return(list(val_int = mean(y), vec_coef = rep(0, ncol(x)), val_threshold = 0))
  }
  
  prev_threshold <- stats::quantile(y, probs = initial_quantile)
  idx <- which(y <= prev_threshold)
  prev_midpoint <- sum(weights[idx]*y[idx])/sum(weights[idx])
  res_glm <- list(val_int = 0, vec_coef = rep(0, ncol(x)))
  iter <- 1
  
  while(iter <= num_iterations){
    # update the regression
    idx <- which(y > prev_threshold)
    if(length(idx) > ncol(x)){
      res_glm <- .glmnet_fancy2(x[idx,,drop = F], y[idx], weights[idx])
    } else {
      return(list(val_int = res_glm$val_int, 
                  vec_coef = res_glm$vec_coef,
                  val_threshold = prev_threshold,
                  val_midpoint = prev_midpoint))
    }
    
    # update the threshold
    tmp <- .update_threshold_glmnet2(x, y, weights, res_glm)
    next_threshold <- tmp$val_threshold
    next_midpoint <- tmp$val_midpoint
    
    if(abs(next_threshold - prev_threshold) <= tol) break()
    iter <- iter+1
  }
  
  list(val_int = res_glm$val_int, vec_coef = res_glm$vec_coef,
       val_threshold = next_threshold)
}

.update_threshold_glmnet2 <- function(x, y, weights, res_glm, tol = 1e-6){
  if(stats::sd(weights*y/sum(weights)) <= 1e-6){
    return(list(val_threshold = max(y), 
                val_midpoint = sum(weights*y)/sum(weights)))
  }
  
  pred_y <- as.numeric(x %*% res_glm$vec_coef) + res_glm$val_int
  
  f <- function(val_threshold, y, pred_y){
    idx <- which(pred_y > val_threshold)
    if(length(idx) == 0) {
      sum(weights*(y - pred_y)^2)
    } else if (length(idx) == length(y)){
      sum(weights*(y - sum(weights * y)/sum(weights))^2)
    } else {
      val1 <- sum(weights[idx]*(y[idx] - pred_y[idx])^2)
      
      y2 <- y[-idx]
      weights2 <- weights[-idx]
      val2 <- sum(weights2*(y2 - sum(weights2*y2)/sum(weights2))^2)
      val1 + val2
    }
  }
  
  val_threshold <- stats::optimize(f, interval = c(min(y), max(y)), maximum = F,
                  y = y, pred_y = pred_y)$minimum
  if(any(y <= val_threshold)){
    idx <- which(pred_y <= val_threshold)
    y2 <- y[idx]
    weights2 <- weights[idx]
    val_midpoint <- sum(weights2*y2)/sum(weights2)
  } else {
    val_midpoint <- 0
  }
  
  list(val_threshold = val_threshold, 
       val_midpoint = val_midpoint)
}

#################################

.transform_est_matrix2 <- function(list_res, ht_map, p1){
  p2 <- length(list_res)

  mat_g <- matrix(0, nrow = p1, ncol = p2)
  vec_g <- rep(0, length = p2)
  vec_threshold <- rep(0, length = p2)
  vec_midpoint <- rep(0, length = p2)
  
  for(j in 1:p2){
    vec_g[j] <- list_res[[j]]$val_int
    vec_threshold[j] <- list_res[[j]]$val_threshold
    vec_midpoint[j] <- list_res[[j]]$val_midpoint
    
    idx_x <- ht_map[[as.character(j)]]
    mat_g[idx_x,j] <- list_res[[j]]$vec_coef
  }
  
  list(mat_g = mat_g, 
       vec_g = vec_g, 
       vec_threshold = vec_threshold,
       vec_midpoint = vec_midpoint)
}

.predict_yfromx2 <- function(mat_x, res_g){
  p2 <- ncol(res_g$mat_g)
  nat_param <- mat_x %*% res_g$mat_g
  
  for(j in 1:p2){
    nat_param[,j] <- nat_param[,j] + res_g$vec_g[j]
  }
  
  if("vec_threshold" %in% names(res_g)){
    for(j in 1:p2){
      threshold <- res_g$vec_threshold[j]
      val <- res_g$vec_midpoint[j]
      res[res[,j] <= threshold, j] <- val
    }
  }
  
  pmax(res, 0)
}
