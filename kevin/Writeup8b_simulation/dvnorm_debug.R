rm(list=ls())
set.seed(10)
p <- 4; n <- 20
mean_vec <- rep(0, p)
sigma <- diag(p)
x_mat <- MASS::mvrnorm(n, mu = mean_vec, Sigma = sigma)
res1 <- .dmvnorm_many_samples(mean = mean_vec,
                              sigma = sigma,
                              x_mat = x_mat)

res2 <- sum(sapply(1:n, function(i){
  .dmvnorm(x = x_mat[i,],
           mean = mean_vec,
           sigma = sigma,
           log = TRUE)
}))

expect_true(abs(res1 - res2) <= 1e-5)

###############

p <- ncol(sigma)
n <- nrow(x_mat)
sigma_inv <- solve(sigma)
determinant_value <- determinant(sigma,
                                 logarithm = TRUE)$modulus
determinant_value <- as.numeric(determinant_value)

x_mat2 <- sweep(x_mat, 
               MARGIN = 2,
               STATS = mean_vec, 
               FUN = "-")
lhs <- x_mat2 %*% sigma_inv
rss <- sum(sapply(1:p, function(j){
  lhs[,j] %*% x_mat2[,j]
}))

alt_sum <- sum(sapply(1:n, function(i){
  t(x_mat[i,]) %*% sigma_inv %*% x_mat[i,]
}))

n * (- 0.5 * determinant_value - 0.5 * p * log(2 * pi)) - 0.5 * rss
