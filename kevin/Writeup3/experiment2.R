rm(list=ls())
k <- 2
n <- 100; p <- 50
cov_mat <- matrix(0, p, p)
cov_mat[1:(p/2), 1:(p/2)] <- 0.5
cov_mat[(p/2+1):p, (p/2+1):p] <- 0.5
diag(cov_mat) <- 1
mat <- MASS::mvrnorm(n, mu = rep(1,p), Sigma = cov_mat)

vec_mean <- matrixStats::colMeans2(mat)
vec_sd <- matrixStats::colSds(mat)
svd_res <- irlba::irlba(mat, nv = k,  
                        center = vec_mean, scale = vec_sd)
zz1 <- svd_res$u %*% diag(svd_res$d/svd_res$d[1])

tmp <- stats::prcomp(mat, center = T, scale. = T)
zz2 <- tmp$x[,1:k]
zz2 <- zz2/svd(zz2)$d[1]

head(zz1)
head(zz2)

########

mat2 <- scale(mat, center = T, scale = T)
tmp <- svd(mat2)
zz3 <- tmp$u[,1:k] %*% diag(tmp$d[1:k]/tmp$d[1])
head(zz3)

zz4 <- mat2 %*% 



