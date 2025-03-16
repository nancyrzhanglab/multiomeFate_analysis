# https://arxiv.org/pdf/1603.05766
.multtest_correction <- function(
    pvalue_vec,
    null_names
){
  stopifnot(length(names(pvalue_vec)) == length(pvalue_vec), 
            all(null_names %in% names(pvalue_vec)))
  
  n <- length(pvalue_vec)
  m <- length(null_names)
  pvalue_vec2 <- rep(NA, n)
  names(pvalue_vec2) <- names(pvalue_vec)
  
  null_values <- pvalue_vec[null_names]
  for(i in 1:n){
    numerator <- length(which(pvalue_vec <= pvalue_vec[i]))+1
    pvalue_vec2[i] <- numerator/(m+1)
  }
  
  return(pvalue_vec2)
}