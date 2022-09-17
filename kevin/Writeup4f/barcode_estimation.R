barcode_estimation <- function(library_vec,
                               lin_mat,
                               dataset_vec,
                               bool_shortcut = T,
                               max_iter = 50,
                               tol = 1e-3,
                               verbose = 1){
  stopifnot(length(library_vec) == length(dataset_vec),
            is.factor(dataset_vec),
            length(library_vec) == ncol(lin_mat),
            inherits(lin_mat, "dgCMatrix"))
  
  if(verbose > 0) print("Initializing multiplier")
  n <- ncol(lin_mat); p <- nrow(lin_mat)
  lin_mat_t <- Matrix::t(lin_mat)
  lineage_count <- Matrix::rowSums(lin_mat)
  
  mean_count_vec <- sapply(levels(dataset_vec), function(dataset){
    idx <- which(dataset_vec == dataset)
    
    mean(Matrix::rowSums(lin_mat[,idx]))
  })
  names(mean_count_vec) <- levels(dataset_vec)
  mean_val <- mean(mean_count_vec)
  multipler_vec <- rep(NA, n)
  for(dataset in levels(dataset_vec)){
    idx <- which(dataset_vec == dataset)
    
    multipler_vec[idx] <- mean_count_vec[dataset]/mean_val
  }
  
  if(verbose > 0) print("Initializing dominant lineages")
  rescaled_mat <- .mult_mat_vec(lin_mat, 1/(multipler_vec*library_vec))
  dominant_assignment_pairs <- lapply(1:n, function(i){
    if(verbose > 1 && i %% floor(n/10) == 0) cat('*')
    
    idx <- .nonzero_col(rescaled_mat, col_idx = i, bool_value = F)
    val <- .nonzero_col(rescaled_mat, col_idx = i, bool_value = T)
    
    if(length(idx) == 0) return(numeric(0))
    if(length(idx) == 1) return(c(idx[1],i))
    if(length(which(val == max(val))) == 1) return(c(idx[which.max(val)],i))
    return(numeric(0))
  })
  dominant_assignment_pairs <- do.call(rbind, dominant_assignment_pairs)
  dominant_mat <- Matrix::sparseMatrix(i = dominant_assignment_pairs[,1],
                                       j = dominant_assignment_pairs[,2],
                                       x = rep(1, nrow(dominant_assignment_pairs)),
                                       dims = c(p,n))
  dominant_mat_t <- Matrix::t(dominant_mat)
  
  if(verbose > 0) print("Initializing betas")
  beta_mat <- sapply(1:p, function(j){
    if(verbose > 1 && j %% floor(p/10) == 0) cat('*')
    
    idx <- .nonzero_col(lin_mat_t, col_idx = j, bool_value = F)
    idx_dominant <- .nonzero_col(dominant_mat_t, col_idx = j, bool_value = F)
    if(length(idx) == 0) return(c(NA, NA))
    
    stopifnot(all(idx_dominant %in% idx))
    if(length(idx_dominant) == 0){
      beta0 <- Matrix::mean(rescaled_mat[j,]) # include even lineages with 0 count
    } else {
      beta0 <- Matrix::mean(rescaled_mat[j,-idx_dominant]) # include even lineages with 0 count
    }
    if(length(idx_dominant) == 0) {
      beta1 <- beta0
    } else {
      beta1 <- max(mean(rescaled_mat[j,idx_dominant]), beta0)
    }
    
    c(beta0, beta1)
  })
  
  beta0 <- beta_mat[1,]; beta0[is.na(beta0)] <- min(beta0, na.rm = T)/2
  beta0[beta0 == 0] <- min(beta0[beta0 != 0])
  beta1 <- beta_mat[2,]; beta1[is.na(beta1)] <- min(beta1, na.rm = T)/2
  pi1 <- 1/p; pi0 <- 1-1/p
  nonzero_list <- sapply(1:n, function(i){
    .nonzero_col(rescaled_mat, col_idx = i, bool_value = F)
  })
  lin_mat <- as.matrix(lin_mat)
  
  iter <- 1
  while(iter <= max_iter){
    if(verbose > 0) print(paste0("On iteration ", iter))
    if(verbose > 0) print("E step")
    posterior_mat1 <- .e_step(beta0 = beta0, 
                              beta1 = beta1,
                              library_vec = library_vec,
                              lin_mat = lin_mat, 
                              multipler_vec = multipler_vec,
                              nonzero_list = nonzero_list,
                              pi0 = pi0,
                              pi1 = pi1,
                              verbose = verbose)
    
    if(verbose > 0) print("M step")
    beta_res <- .m_step(beta0 = beta0, 
                        library_vec = library_vec,
                        lin_mat = lin_mat, 
                        multipler_vec = multipler_vec,
                        posterior_mat1 = posterior_mat1,
                        bool_shortcut = bool_shortcut,
                        verbose = verbose)
    beta0 <- beta_res$beta0; beta1 <- beta_res$beta1
    iter <- iter + 1
  }
  
  list(beta0 = beta0,
       beta1 = beta1,
       dominant_mat = dominant_mat,
       posterior_mat1 = posterior_mat1)
}

##################

.e_step <- function(beta0, 
                    beta1,
                    library_vec,
                    lin_mat, 
                    multipler_vec,
                    nonzero_list,
                    pi0,
                    pi1,
                    verbose = 1){
  stopifnot(is.matrix(lin_mat))
  
  n <- ncol(lin_mat)
  posterior_mat1 <- sapply(1:n, function(i){
    if(verbose > 1 && i %% floor(n/10) == 0) cat('*')
    
    vec0 <- stats::dpois(x = lin_mat[,i],
                         lambda = beta0*library_vec[i]*multipler_vec[i])
    vec1 <- stats::dpois(x = lin_mat[,i],
                         lambda = beta1*library_vec[i]*multipler_vec[i])
    posterior_vec1 <- vec1*pi1/(vec1*pi1 + vec0*pi0)
    
    if(length(nonzero_list[[i]]) > 0){
      posterior_vec1[-nonzero_list[[i]]] <- pmin(posterior_vec1[-nonzero_list[[i]]], 0.5)
    } else {
      posterior_vec1 <- pmin(posterior_vec1, 0.5)
    }
    
    posterior_vec1
  })
  
  posterior_mat1
}

.m_step <- function(beta0, 
                    library_vec,
                    lin_mat, 
                    multipler_vec,
                    posterior_mat1,
                    bool_shortcut,
                    tol = 1e-4,
                    verbose = 1){
  stopifnot(is.matrix(lin_mat))
  
  p <- nrow(lin_mat)
  beta_mat <- sapply(1:p, function(j){
    if(verbose >= 3) print(paste0(j , " of ", p))
    if(verbose == 2 && j %% floor(p/10) == 0) cat('*')
    
    nonzero_idx <- which(posterior_mat1[j,] < 1-tol)
    if(length(nonzero_idx) == 0){
      beta0_val <- beta0[j]
    } else {
      df <- data.frame(lib = lin_mat[j,nonzero_idx], 
                       covariate = multipler_vec[nonzero_idx] * library_vec[nonzero_idx])
      
      if(bool_shortcut){
        weighted_num <- sum((1-posterior_mat1[j,nonzero_idx])*(log1p(df[,"lib"])/df[,"covariate"]))
        weighted_denom <- sum((1-posterior_mat1[j,nonzero_idx]))
        beta0_val <- exp(weighted_num/weighted_denom)-1
      } else {
        glm_res0 <- suppressWarnings(stats::glm(lib ~ . - 1, data = df, family = stats::poisson, 
                                                weights = 1-posterior_mat1[j,nonzero_idx]))
        beta0_val <- stats::coef(glm_res0)[1]
      }
      if(beta0_val <= 0) beta0_val <- beta0[j]
    }
    
    nonzero_idx <- which(posterior_mat1[j,] > tol)
    if(length(nonzero_idx) == 0){
      beta1_val <- beta0_val
    } else {
      df <- data.frame(lib = lin_mat[j,nonzero_idx], 
                       covariate = multipler_vec[nonzero_idx] * library_vec[nonzero_idx])
      if(bool_shortcut){
        weighted_num <- sum(posterior_mat1[j,nonzero_idx]*(log1p(df[,"lib"])/df[,"covariate"]))
        weighted_denom <- sum(posterior_mat1[j,nonzero_idx])
        beta1_val <- exp(weighted_num/weighted_denom)-1
      } else {
        glm_res1 <- suppressWarnings(stats::glm(lib ~ . - 1, data = df, family = stats::poisson, 
                                                weights = posterior_mat1[j,nonzero_idx]))
        beta1_val <- stats::coef(glm_res0)[1]
      }
      if(beta1_val <= 0) beta1_val <- beta0_val
    }
    
    c(beta0_val, beta1_val)
  })
  
  df <- as.data.frame(t(log(beta_mat))); colnames(df) <- c("beta0", "beta1")
  lm_res <- stats::lm(beta0 ~ beta1, data = df)
  beta1 <- beta_mat[2,]
  beta0 <- exp(lm_res$fitted.values)
  beta0 <- pmax(pmin(beta0, beta1), min(beta_mat[1,]))
  stopifnot(all(beta1 > 0), all(beta0 > 0), all(beta0 <= beta1))
  
  list(beta0 = beta0,
       beta1 = beta1)
}

##################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}

