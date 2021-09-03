# see https://github.com/ShobiStassen/VIA/blob/master/VIA/core.py#L1094
# https://igraph.org/r/doc/page_rank.html
form_transition <- function(snn, 
                            lazy_param = 0.85,
                            teleport_param = 0.99){
  n <- nrow(snn)
  P <- snn
  
  for(i in 1:n){
    idx <- which(P[i,] != 0)
    min_val <- min(P[i,idx])
    max_val <- max(P[i,idx])
    P[i,idx] <- exp(-(P[i,idx]- min_val)/max_val)
  }
  P <- diag(1/matrixStats::rowSums2(P)) %*% P
  
  P <- lazy_param*P + (1-lazy_param)*diag(ncol(P))
  P <- teleport_param*P + (1-teleport_param)*matrix(1/n, n, n)
  
  rownames(P) <- rownames(adj_mat)
  colnames(P) <- colnames(adj_mat)
  
  P
}

# try the diffusion distance, see 
# https://www.stat.berkeley.edu/~mmahoney/s15-stat260-cs294/Lectures/lecture15-12mar15.pdf
# https://academic.oup.com/bioinformatics/article/31/18/2989/241305
# https://arxiv.org/pdf/1406.0013.pdf
# https://arxiv.org/pdf/0811.0121.pdf
# https://mathworld.wolfram.com/LeftEigenvector.html
# https://mathworld.wolfram.com/Eigenvector.html
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
    stopifnot(sum(abs(sort(right_eigen$values) - sort(left_eigen$values))) <= 1e-6)
    
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
    sum((c(eigenvalues[-1])^x*(right_vector[idx1,-1] - right_vector[idx2,-1]))^2)
  })))
}