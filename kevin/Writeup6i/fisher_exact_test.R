# https://github.com/linnykos/permanent_notes/blob/master/test_classical/Fishers-Exact-test-for-two-Proportions.pdf
exact_fisher_pvalue <- function(count_group1,
                                count_group2){
  stopifnot(length(count_group1) == 2,
            length(count_group2) == 2)
  if(length(names(count_group1)) > 0){
    stopifnot(all(names(count_group1) == names(count_group2)))
  }
  
  n1 <- sum(count_group1)
  n2 <- sum(count_group2)
  x1 <- count_group1[1]
  x2 <- count_group2[1]
  N <- n1+n2
  m <- x1+x2
  
  test_stat <- .compute_teststat(m, n1, n2, N, x1, x2)
  
  splits <- .enumerate_splits(m)
  other_values <- sapply(1:nrow(splits), function(i){
    .compute_teststat(m, n1, n2, N, 
                      x1 = splits[i,"x1"], 
                      x2 = splits[i,"x2"])
  })
  bool_idx <- which(other_values >= test_stat)
  logprob_val <- .add_logged_probs(flip_sign = T, 
                                   log_vec = other_values[bool_idx])
  log10prob_val <- logprob_val/(log(10))
  log10prob_val <- -log10prob_val
  
  list(log10prob_val = log10prob_val,
       test_stat = test_stat)
}

.compute_teststat <- function(m, n1, n2, N, x1, x2){
  if(x1 > n1 | x2 > n2) return(NA)
  -(lchoose(n1,x1) + lchoose(n2,x2) - lchoose(N,m))
}

.enumerate_splits <- function(m){
  mat <- cbind(0:m, m:0)
  colnames(mat) <- c("x1", "x2")
  mat
}

# compute log(exp(log_val1) + exp(log_val2) + ...)
.add_logged_probs <- function(flip_sign, log_vec){
  if(flip_sign) log_vec <- -log_vec
  c <- min(log_vec)
  log(sum(exp(log_vec + c))) - c
}

#################

# for debugging purposes only
exact_fisher_pvalue_simple <- function(count_group1,
                                       count_group2){
  stopifnot(length(count_group1) == 2,
            length(count_group2) == 2)
  if(length(names(count_group1)) > 0){
    stopifnot(all(names(count_group1) == names(count_group2)))
  }
  
  n1 <- sum(count_group1)
  n2 <- sum(count_group2)
  x1 <- count_group1[1]
  x2 <- count_group2[1]
  N <- n1+n2
  m <- x1+x2
  
  test_stat <- .compute_teststat2(m, n1, n2, N, x1, x2)
  
  splits <- .enumerate_splits(m)
  other_values <- sapply(1:nrow(splits), function(i){
    .compute_teststat2(m, n1, n2, N, 
                      x1 = splits[i,"x1"], 
                      x2 = splits[i,"x2"])
  })
  bool_idx <- which(other_values <= test_stat)
  prob_val <- sum(other_values[bool_idx])
  
  list(other_values = other_values,
       prob_val = prob_val,
       test_stat = -log(test_stat))
}

.compute_teststat2 <- function(m, n1, n2, N, x1, x2){
  if(x1 > n1 | x2 > n2) return(NA)
  choose(n1,x1)*choose(n2,x2)/choose(N,m)
}
