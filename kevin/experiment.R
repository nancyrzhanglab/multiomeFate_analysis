tmp <- theta[j,]
V <- theta[combn_mat[,i],]

n <- nrow(V); K <- ncol(V)

idx <- n
theta2 <- tmp - V[idx,]; V2 <- V[-idx,,drop = F] - rep(V[idx,], each = n-1)

D <- tcrossprod(V2)
d <- V2 %*% tmp
b0 <- c(rep(0,n-1), -1)
A <- cbind(diag(n-1), rep(-1, n-1))

obj <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b0)
value <- sqrt(max(2*obj$value + .l2norm(tmp)^2, 0))

comb <- rep(0, n)
comb[-idx] <- obj$solution; comb[idx] <- 1-sum(obj$solution)
