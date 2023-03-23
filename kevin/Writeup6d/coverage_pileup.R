coverage_pileup <- function(
    object,
    gene, # name of gene
    assay = "ATAC",
    cells = NULL, #NULL or names of cells
    extend.downstream = 5000,
    extend.upstream = 5000,
    sep = c("-", "-"),
    window = 100
){
  cutmat <- .compute_cutmatrix(
    object = object,
    gene = gene,
    assay = assay,
    cells = cells, 
    extend.downstream = extend.downstream,
    extend.upstream = extend.upstream,
    sep = sep,
    window = window
  )
  if(all(is.null(cutmat))) return(NULL)
  
  # find all the peak regions
  peak_matrix <- extract_peaks(
    object = object,
    gene = gene, 
    assay = assay,
    extend.downstream = extend.downstream,
    extend.upstream = extend.upstream,
    sep = sep
  )
  
  # compute the column indices of the midpoints
  peak_midpoint <- apply(peak_matrix, 1, function(x){round(mean(x))})
  
  # for each fragment (in any of the cells), compute the signed distance to the nearest peak midpoint
  # - WARNING: Technically, the start-end of each fragment might be split into two different nearest midpoints...
  res <- lapply(1:ncol(cutmat), function(j){
    val <- length(.nonzero_col(mat = cutmat,
                        col_idx = j,
                        bool_value = T))
    if(val == 0) return(numeric(0))
    
    current_bp <- as.numeric(colnames(cutmat)[j])
    idx <- which.min(abs(current_bp - peak_midpoint))
    distance <- current_bp - peak_midpoint[idx]
    c(val, distance)
  })
  names(res) <- colnames(cutmat)
  res <- res[which(sapply(res, length) > 0)]
  
  # return 2-column vector
  # - one vector is the number of fragments in a particular basepair point
  # - one vector of is distances to the nearest peak region
  res2 <- do.call(res, rbind)
  rownames(res2) <- names(res)
  colnames(res2) <- c("value", "distance")
  res2
  
}

#######################

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
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
