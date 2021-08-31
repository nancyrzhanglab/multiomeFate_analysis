form_transition <- function(adj_mat, normalize = T){
  P <- adj_mat
  tmp <- matrixStats::rowSums2(P)
  if(normalize) P <- diag(1/tmp) %*% P %*% diag(1/tmp)
  diag(1/matrixStats::rowSums2(P)) %*% P
}

extract_eigen <- function(P, check = F){
  right_eigen <- eigen(P, symmetric = F)
  left_eigen <- eigen(t(P), symmetric = F)
  if(any(left_eigen$vectors[,1] < 0)) left_eigen$vectors[,1] <- -left_eigen$vectors[,1]
  eigenvalues <- left_eigen$values 
  
  if(check){
    n <- nrow(P)
    stopifnot(sum(abs(P%*%right_eigen$vectors[,2] - right_eigen$values[2]*right_eigen$vectors[,2])) <= 1e-6)
    stopifnot(sum(abs(left_eigen$vectors[,2]%*%P - left_eigen$values[2]*left_eigen$vectors[,2])) <= 1e-6)
    stopifnot(sum(abs(right_eigen$vectors[,1]+rep(1/sqrt(n)))) <= 1e-6 |
                sum(abs(right_eigen$vectors[,1]-rep(1/sqrt(n)))) <= 1e-6)
    stopifnot(sum(abs(right_eigen$values - left_eigen$values)) <= 1e-6)
    
    stopifnot(abs(sum(left_eigen$vectors[,2]^2) - 1) <= 1e-6)
  }
  
  # normalize
  left_vector <- left_eigen$vectors
  for(i in 2:ncol(left_vector)){
    left_vector[,i] <- left_vector[,i]*sqrt(left_vector[,1])
  }
  right_vector <- right_eigen$vectors
  for(i in 1:ncol(right_vector)){
    right_vector[,i] <- right_vector[,i]/sqrt(left_vector[,1])
  }
  
  list(eigenvalues = eigenvalues, 
       left_vector = left_vector,
       right_vector = right_vector)
}

diffusion_distance <- function(eigenvalues, right_vector, 
                               idx1, idx2, 
                               time_vec = 1:min(40, length(eigenvalues))){
  sqrt(sum(sapply(time_vec, function(x){
    eigenvalues^(2*x)*sum((right_vector[idx1,-1] - right_vector[idx2,-1])^2)
  })))
}