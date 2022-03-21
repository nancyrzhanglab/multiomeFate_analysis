
.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}


.form_snn_mat <- function(mat, 
                          num_neigh,
                          bool_cosine, # suggested: TRUE
                          bool_intersect, # suggested: TRUE
                          tol = 1e-4,
                          verbose = 0){
  stopifnot(num_neigh >= min_deg, min_deg >= 0)
  
  if(bool_cosine) {
    l2_vec <- apply(mat, 1, .l2norm)
    l2_vec[l2_vec <= tol] <- tol
    mat <- .mult_vec_mat(1/l2_vec, mat)
  }
  
  if(verbose >= 3) print("Compute NNs")
  n <- nrow(mat)
  nn_mat <- RANN::nn2(mat, k = num_neigh+1)$nn.idx
  if(all(nn_mat[,1] == 1:n)){
    nn_mat <- nn_mat[,-1,drop = F]
  }
  
  if(verbose >= 3) print("Forming NN matrix")
  i_vec <- rep(1:n, times = ncol(nn_mat))
  j_vec <- as.numeric(nn_mat)
  
  sparse_mat <- Matrix::sparseMatrix(i = i_vec,
                                     j = j_vec,
                                     x = rep(1, length(i_vec)),
                                     dims = c(n,n),
                                     repr = "C")
  
  if(bool_intersect) {
    sparse_mat <- sparse_mat * Matrix::t(sparse_mat)
  } else {
    sparse_mat <- sparse_mat + Matrix::t(sparse_mat)
    sparse_mat@x <- rep(1, length(sparse_mat@x))
  }
  
  sparse_mat
}

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