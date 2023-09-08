trials <- 1000
n_vec <- round(seq(2, 100, length.out = 20))
mat <- matrix(NA, nrow = length(n_vec), ncol = trials)
rownames(mat) <- paste0("n:", n_vec)
  
set.seed(10)
for(i in 1:length(n_vec)){
  mat[i,] <- sapply(1:trials, function(trial){
    n <- n_vec[i]
    vec <- stats::rpois(n, lambda = 1)
    (sum(vec) - n)/n
  })
}

df <- cbind(as.numeric(mat), rep(n_vec, times = trials))
plot(df[,2], df[,1], pch = 16, col = rgb(0.5,0.5,0.5,0.2))

########################################

mat2 <- matrix(NA, nrow = length(n_vec), ncol = trials)
rownames(mat2) <- paste0("n:", n_vec)

set.seed(10)
for(i in 1:length(n_vec)){
  mat2[i,] <- sapply(1:trials, function(trial){
    n <- n_vec[i]
    vec <- stats::rpois(n, lambda = 1)
    (sum(vec) - n)/sqrt(n) #sqrt(n) is the standard deviaion
  })
}

df2 <- cbind(as.numeric(mat2), rep(n_vec, times = trials))
plot(df2[,2], df2[,1], pch = 16, col = rgb(0.5,0.5,0.5,0.2))
