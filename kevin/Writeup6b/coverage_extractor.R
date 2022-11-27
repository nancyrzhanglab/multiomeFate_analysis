coverage_extractor <- function(
    object,
    region,
    which_ident,
    assay = "ATAC",
    cells = NULL,
    extend.downstream = 1000,
    extend.upstream = 1000,
    group.by = NULL,
    idents = NULL,
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
  if (!is.null(x = idents)) {
    ident.cells <- Signac:::WhichCells(object = object, idents = idents)
    cells <- intersect(x = cells, y = ident.cells)
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
    idents = idents
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
  group.scale.factors <- group.scale.factors[[1]]
  scale.factor <- Signac::SetIfNull(
    x = scale.factor, y = median(x = group.scale.factors)
  )
  scale.factor <- scale.factor[[1]]
 
  ##########
  
  levels.use <- levels(x = obj.groups)
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
  
  # multiply cutmat by the rolling average matrix
  avg_mat <- .construct_averaging_mat(p = ncol(cutmat),
                                      window = window)
  cutmat_list <- lapply(cutmat_list, function(mat){
    tmp <- mat %*% avg_mat
    colnames(tmp) <- colnames(mat)[floor(window/2):(ncol(mat)-floor(window/2))]
    tmp
  })
  
  # normalize all the value
  for(i in 1:length(cutmat_list)){
    cutmat_list[[i]]@x <- cutmat_list[[i]]@x / group.scale.factors[which_ident[i]] * scale.factor
  }
  
  # compute the median and the quantiles
  coverage_mean <- sapply(cutmat_list, function(mat){
    Matrix::colMeans(mat)
  })
  
  # vv=== old code ===vv #
  coverages <- Signac:::ApplyMatrixByGroup(
    mat = cutmat,
    fun = colSums,
    groups = obj.groups,
    group.scale.factors = group.scale.factors,
    scale.factor = scale.factor,
    normalize = TRUE
  )
  coverages <- dplyr::group_by(.data = coverages, group)
  coverages <- dplyr::mutate(.data = coverages, 
                             coverage = RcppRoll::roll_sum(
                               x = norm.value, 
                               n = window, 
                               fill = NA, 
                               align = "center"
                             ))
  coverages <- dplyr::ungroup(x = coverages)
  coverages <- coverages[!is.na(x = coverages$coverage), ]
  coverages <- dplyr::group_by(.data = coverages, group)
  
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