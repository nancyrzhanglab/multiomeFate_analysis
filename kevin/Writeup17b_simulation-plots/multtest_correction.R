.multtest_correction <- function(
    pvalue_vec
){
  stopifnot(length(names(pvalue_vec)) == length(pvalue_vec))
  
  n <- length(pvalue_vec)
  
  z_vec <- Rmpfr::qnormI(pvalue_vec)
  neg_idx <- intersect(which(is.infinite(z_vec)), which(sign(z_vec) < 0))
  if(length(neg_idx) > 0){
    min_val <- stats::quantile(z_vec[-neg_idx], probs = 0.01)
    z_vec <- pmax(z_vec, min_val)
  }
  pos_idx <- intersect(which(is.infinite(z_vec)), which(sign(z_vec) > 0))
  if(length(pos_idx) > 0){
    max_val <- stats::quantile(z_vec[-pos_idx], probs = 0.99)
    z_vec <- pmin(z_vec, max_val)
  }
 
  median_val <- stats::median(z_vec)
  sd_val <- stats::sd(z_vec)
  
  pvalue_vec2 <- sapply(z_vec, function(z){
    Rmpfr::pnorm(z, mean = median_val, sd = sd_val)
  })
  names(pvalue_vec2) <- names(pvalue_vec)
  
  return(pvalue_vec2)
}
