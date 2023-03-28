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
  if(length(res) == 0){
    res2 <- matrix(NA, nrow = 0, ncol = 2)
    colnames(res2) <- c("value", "distance")
  } else {
    # return 2-column vector
    # - one vector is the number of fragments in a particular basepair point
    # - one vector of is distances to the nearest peak region
    res2 <- do.call(rbind, res)
    rownames(res2) <- names(res)
    colnames(res2) <- c("value", "distance")
  }
  
  peak_width <- apply(peak_matrix, 1, diff)
  list(pileup_mat = res2,
       peak_width_median = stats::median(peak_width),
       peak_width_max = max(peak_width))
}

compute_pileup_curve <- function(
    pileup_mat,
    inflation_factor = 2,
    peak_width_max = NULL, # to normalize, peak_width_max should be not NULL
    peak_width_median = NULL,
    window = 100
){
  if(nrow(pileup_mat) == 0){
    cutvec <- rep(0, 201)
    names(cutvec) <- seq(-100,100)
  } else {
    cutvec <- .form_cutvec_from_pileup(pileup_mat,
                                       window = 100)
  }
  sum_mat <- .construct_sum_mat_triangle(p = length(cutvec),
                                         window = window)
  vec <- cutvec %*% sum_mat
  vec <- as.numeric(vec)
  names(vec) <- names(cutvec)[1:length(vec)]
  
  # center the vector
  left_bp <- as.numeric(names(vec)[1])
  right_bp <- as.numeric(names(vec)[length(vec)])
  if(left_bp < 0 & right_bp > 0){
    if(abs(left_bp) != abs(right_bp)){
      if(abs(left_bp) > abs(right_bp)){
        vec <- c(vec, rep(0, abs(left_bp)-abs(right_bp)))
        names(vec) <- seq(left_bp, abs(left_bp))
      } else {
        vec <- c(rep(0, abs(right_bp)-abs(left_bp)), vec)
        names(vec) <- seq(-right_bp, right_bp)
      }
    }
  } else if(left_bp < 0 & right_bp < 0){
    big_num <- abs(left_bp)
    small_num <- abs(right_bp)
    stopifnot(big_num > small_num)
    # print("Scenario A")
    # print(paste0(big_num, ":", small_num, "::", left_bp, ":", right_bp, "::", length(vec)))
    
    vec <- c(vec, rep(0, small_num+big_num))
    # print(length(vec))
    names(vec) <- seq(-big_num, big_num)
  } else if(left_bp > 0 & right_bp > 0) {
    big_num <- abs(right_bp)
    small_num <- abs(left_bp)
    stopifnot(big_num > small_num)
    # print("Scenario B")
    # print(paste0(big_num, ":", small_num, "::", left_bp, right_bp, "::", length(vec)))
    
    vec <- c(rep(0, small_num+big_num), vec)
    # print(length(vec))
    names(vec) <- seq(-big_num, big_num)
  }
 
  # normalize
  if(!is.null(peak_width_max)){
    # figure out which entries are in the "middle"
    midpoint_idx <- round(length(vec)/2)
    midpoint_window <- round(c(-1,1)*inflation_factor*peak_width_max + midpoint_idx)
    midpoint_window <- pmax(pmin(midpoint_window, length(vec)), 1)
    
    # normalize by dividing by the mean value outside the peak region
    idx <- c(1:length(vec))[-(midpoint_window[1]:midpoint_window[2])]
    if(length(idx) > 0) vec <- vec/mean(vec[idx])
  }
  
  # compute score
  if(!is.null(peak_width_median)){
    # figure out which entries are in the "middle"
    midpoint_idx <- round(length(vec)/2)
    midpoint_window <- round(c(-.5,.5)*peak_width_median + midpoint_idx)
    midpoint_window <- pmax(pmin(midpoint_window, length(vec)), 1)
    
    # mean value inside the peak region
    score <- mean(vec[midpoint_window[1]:midpoint_window[2]])
  } else {
    score <- NULL
  }
  
  list(inflation_factor = inflation_factor,
       peak_width_max = peak_width_max,
       peak_width_median = peak_width_median,
       pileup_vec = vec,
       score = score)
}

#######################

.form_cutvec_from_pileup <- function(pileup_mat,
                                     window){
  width <- diff(range(pileup_mat[,"distance"]))
  width <- width + 2*window
  
  min_distance <- min(pileup_mat[,"distance"])
  column_vec <- pileup_mat[,"distance"] - min_distance + 1 + window
  vec <- rep(0, width)
  vec[column_vec] <- pileup_mat[,"value"]
  min_basepair <- min(pileup_mat[,"distance"]) - window
  names(vec) <- seq(min_basepair,(min_basepair+length(vec)-1))
  vec
}

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
