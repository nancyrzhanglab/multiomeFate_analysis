coverage_extractor_singlecell <- function(
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
  
  # normalize values for each cell based on that row
  cell_size <- pmax(Matrix::rowSums(cutmat),1)
  cutmat <- .mult_vec_mat(1/cell_size, cutmat)
  
  # smooth values 
  sum_mat <- .construct_sum_mat_triangle(p = ncol(cutmat),
                                         window = window)
  tmp <- cutmat %*% sum_mat
  colnames(tmp) <- colnames(cutmat)[floor(window/2):(floor(window/2)+ncol(tmp)-1)]
  rownames(tmp) <- cells
  tmp
}

extract_peaks <- function(
    object,
    gene, # name of gene
    assay = "ATAC",
    extend.downstream = 5000,
    extend.upstream = 5000,
    sep = c("-", "-")
){
  
  tmp <- Signac::LookupGeneCoords(
    object = object,
    gene = gene,
    assay = assay
  )
  # make sure gene exists
  if(all(is.null(tmp))) return(NULL)
  
  region <- Signac:::FindRegion(
    object = object,
    region = gene,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  all_atac_peak <- object[["ATAC"]]@ranges
  overlap_res <- GenomicRanges::findOverlaps(
    query = all_atac_peak,
    subject = region
  )
  
  if(length(overlap_res) > 0){
    region_gene_peaks <- all_atac_peak[overlap_res@from]
    for(i in 1:length(region_gene_peaks)){
      region_gene_peaks[i] <- intersect(x = region_gene_peaks[i],
                                        y = region)
    }
  }
  
  ranges_obj <- region_gene_peaks@ranges
  tmp <- cbind(ranges_obj@start, ranges_obj@start + ranges_obj@width - 1)
  colnames(tmp) <- c("start", "end")
  tmp
}

compute_windowsize <- function(
    object,
    gene,
    assay = "ATAC",
    base_length = 1000,
    extend.downstream = 5000,
    extend.upstream = 5000,
    sep = c("-", "-"),
    width_threshold = 100000 #1e5
){
  tmp <- Signac::LookupGeneCoords(
    object = object,
    gene = gene,
    assay = assay
  )
  # make sure gene exists
  if(all(is.null(tmp))) return(NULL)
  
  region <- Signac:::FindRegion(
    object = object,
    region = gene,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  width <- IRanges::end(x = region) - IRanges::start(x = region)
  
  if(width <= width_threshold){ #1e5
    base_length
  } else {
    round(base_length * width/width_threshold)
  }
}

######################################################

.compute_cutmatrix <- function(
    object,
    gene, # name of gene
    assay = "ATAC",
    cells = NULL, #NULL or names of cells
    extend.downstream = 5000,
    extend.upstream = 5000,
    sep = c("-", "-"),
    window = 100
){
  
  cells <- Signac:::SetIfNull(x = cells, y = colnames(x = object))
  
  tmp <- Signac::LookupGeneCoords(
    object = object,
    gene = gene,
    assay = assay
  )
  # make sure gene exists
  if(all(is.null(tmp))) return(NULL)
  
  region <- Signac:::FindRegion(
    object = object,
    region = gene,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  # form cutmatrix
  cutmat <- Signac:::CutMatrix(
    object = object,
    region = region,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  colnames(cutmat) <- (IRanges::start(x = region)):(IRanges::end(x = region))
  
  cutmat
}

.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

.construct_sum_mat_triangle <- function(p,
                                        window){
  stopifnot(window + 1 < p)
  vec <- c(seq(0, 1, length.out = ceiling(window/2)+1),
           seq(1, 0, length.out = ceiling(window/2)+1)[-1])
  vec <- vec[1:window]
  
  i_vec <- unlist(lapply(1:(p-window+1), function(i){
    i:(i+window-1)
  }))
  j_vec <- unlist(lapply(1:(p-window+1), function(j){
    rep(j, window)
  }))
  x_vec <- unlist(lapply(1:(p-window+1), function(x){
    vec
  }))
  
  Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(p, p-window+1))
}
