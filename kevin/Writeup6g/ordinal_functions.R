.compute_confusion_score <- function(prediction_vec, truth_vec){
  uniq_val <- sort(unique(truth_vec))
  stopifnot(all(prediction_vec %in% uniq_val), 
            length(prediction_vec) == length(truth_vec))
  
  k <- length(uniq_val)
  len_vec <- sapply(uniq_val, function(val){
    length(which(truth_vec == val))
  })
  
  tmp <- table(prediction_vec, truth_vec)
  tmp <- tmp[,uniq_val]
  conf_mat <- matrix(0, nrow = k, ncol = k)
  rownames(conf_mat) <- uniq_val
  colnames(conf_mat) <- uniq_val
  conf_mat[rownames(tmp),] <- tmp
  stopifnot(nrow(conf_mat) == ncol(conf_mat))
  
  conf_mat <- conf_mat %*% diag(1/len_vec)
  mean(diag(conf_mat))
}

.five_fold_cv <- function(x_mat, y_vec){
  stopifnot(nrow(x_mat) == length(y_vec),
            is.factor(y_vec), all(diff(as.numeric(levels(y_vec))) > 0))
  
  n <- nrow(x_mat)
  uniq_val <- sort(unique(y_vec))
  k <- length(uniq_val)
  len_vec <- sapply(uniq_val, function(val){
    length(which(y_vec == val))
  })
  names(len_vec) <- uniq_val
  
  # construct weights
  idx_list <- lapply(uniq_val, function(val){
    which(y_vec == val)
  })
  weight_vec <- rep(NA, n)
  for(i in 1:k){
    weight_vec[idx_list[[i]]] <- 1/length(idx_list[[i]])
  }
  
  # construct folds
  fold_vec <- rep(NA, n)
  for(i in 1:k){
    tmp_vec <- rep(1:5, each = ceiling(length(idx_list[[i]])/5))
    tmp_vec <- sample(tmp_vec)[1:length(idx_list[[i]])]
    fold_vec[idx_list[[i]]] <- tmp_vec
  }
  
  # apply cross-validation
  cv_score <- sapply(1:5, function(fold){
    fold_idx <- which(fold_vec == fold)
    train_x <- x_mat[-fold_idx,,drop = F]
    test_x <- x_mat[fold_idx,,drop = F]
    train_y <- y_vec[-fold_idx]
    test_y <- y_vec[fold_idx]
    weight_tmp <- weight_vec[-fold_idx]
    
    train_df <- cbind(train_y, train_x)
    colnames(train_df)[1] <- "y"
    train_df <- as.data.frame(train_df)
    train_df[,"y"] <- as.factor(train_df[,"y"])
    
    oridinal_res <- ordinal::clm(y ~ ., data = train_df, weights = weight_tmp)
    pred_vec <- stats::predict(oridinal_res, newdata = as.data.frame(test_x))$fit
    pred_vec <- sapply(1:nrow(pred_vec), function(i){which.max(pred_vec[i,])})
    
    .compute_confusion_score(
      prediction_vec = pred_vec,
      truth_vec = test_y
    )
  })
  
  mean(cv_score)
}

.permutation_null_score <- function(x_mat, y_vec, 
                                    trials = 100,
                                    verbose = 1){
  n <- nrow(x_mat)
  
  cv_score_trials <- sapply(1:trials, function(trial){
    if(verbose > 0 && trial %% floor(trials/10) == 0) cat('*')
    y_vec_null <- sample(y_vec)
    .five_fold_cv(x_mat = x_mat,
                  y_vec = y_vec_null)
  })
  
  cv_score_trials
}