coverage_extractor <- function(
    object,
    region,
    group.by, # name of ident you want to use
    which_ident, # which values of the ident do you want to output cells for?
    assay = "ATAC",
    cells = NULL,
    extend.downstream = 1000,
    extend.upstream = 1000,
    multipler = 1e6,
    scale.factor = NULL,
    sep = c("-", "-"),
    window = 100
){
  
  # some default things
  cells <- Signac:::SetIfNull(x = cells, y = colnames(x = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("Requested assay is not a ChromatinAssay.")
  }
  if (!is.null(x = group.by)) {
    Seurat::Idents(object = object) <- group.by
  }

  ##########
  
  region <- Signac:::FindRegion(
    object = object,
    region = region,
    sep = sep,
    assay = assay[[1]],
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  cells.per.group <- Signac:::CellsPerGroup(
    object = object,
    group.by = group.by
  )
  obj.groups <- Signac:::GetGroups(
    object = object,
    group.by = group.by,
    idents = NULL
  )
  cm.list <- list()
  sf.list <- list()
  gsf.list <- list()
  
  reads.per.group <- Signac:::AverageCounts(
    object = object,
    assay = assay,
    group.by = group.by,
    verbose = FALSE
  )
  cutmat <- Signac:::CutMatrix(
    object = object,
    region = region,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  colnames(cutmat) <- (IRanges::start(x = region)):(IRanges::end(x = region))
  group.scale.factors <- suppressWarnings(reads.per.group * cells.per.group)
  scale.factor <- Signac:::SetIfNull(
    x = scale.factor, y = median(x = group.scale.factors)
  )
  
  ##########
  
  chromosome <- as.character(x = seqnames(x = region))
  start.pos <- IRanges::start(x = region)
  end.pos <- IRanges::end(x = region)
  
  # seperate cutmat by the relevant rows
  idx_list <- lapply(which_ident, function(ident){
    which(Seurat::Idents(object) == ident)
  })
  cutmat_list <- lapply(idx_list, function(idx_vec){
    cutmat[idx_vec,]
  })
  
  ##############################
  # compute the normalized cut matrix
  ##############################
  
  # multiply cutmat by the rolling average matrix to smooth out the entries
  avg_mat <- .construct_averaging_mat(p = ncol(cutmat),
                                      window = window)
  cutmat_norm_list <- lapply(cutmat_list, function(mat){
    tmp <- mat %*% avg_mat
    colnames(tmp) <- colnames(mat)[floor(window/2):(ncol(mat)-floor(window/2))]
    tmp
  })
  
  # normalize all the value
  for(i in 1:length(cutmat_norm_list)){
    cutmat_norm_list[[i]]@x <- cutmat_norm_list[[i]]@x / group.scale.factors[which_ident[i]] * scale.factor
  }
  
  # compute the median and the quantiles
  coverage_mean <- sapply(cutmat_norm_list, function(mat){
    Matrix::colMeans(mat)
  })
  colnames(coverage_mean) <- which_ident
  rownames(coverage_mean) <- colnames(cutmat_norm_list[[1]])
  
  ##############################
  # now, for the unprocessed version, compute the list where you see how many cells have a fragment at a bp
  ##############################
  
  coverage_count <- sapply(cutmat_list, function(mat){
    mat@x <- rep(1, length(mat@x))
    Matrix::colSums(mat)
  })
  colnames(coverage_count) <- which_ident
  coverage_count <- coverage_count[rownames(coverage_mean),]
  
  total_vec <- sapply(cutmat_list, nrow)
  names(total_vec) <- which_ident
  
  list(coverage_mean = coverage_mean,
       coverage_count = coverage_count,
       total_vec = total_vec)
}

.construct_averaging_mat <- function(p,
                                     window){
  stopifnot(window + 1 < p)
  
  i_vec <- unlist(lapply(1:(p-window+1), function(i){
    i:(i+window-1)
  }))
  j_vec <- unlist(lapply(1:(p-window+1), function(j){
    rep(j, window)
  }))
  x_vec <- unlist(lapply(1:(p-window+1), function(x){
    rep(1/window, window)
  }))
  
  Matrix::sparseMatrix(i = i_vec, j = j_vec, x = x_vec, dims = c(p, p-window+1))
}