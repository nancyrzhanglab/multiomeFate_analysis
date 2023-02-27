coverage_extractor <- function(
    object,
    gene,
    assay = "ATAC",
    extend.downstream = 1000,
    extend.upstream = 1000,
    sep = c("-", "-")
){
  
  # some default things
  cells <- Signac:::SetIfNull(x = NULL, y = colnames(x = object))
  if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
    stop("Requested assay is not a ChromatinAssay.")
  }
  
  ##########
  
  region <- Signac:::FindRegion(
    object = object,
    region = gene,
    sep = sep,
    assay = assay[[1]],
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  cutmat <- Signac:::CutMatrix(
    object = object,
    region = region,
    assay = assay,
    cells = cells,
    verbose = FALSE
  )
  colnames(cutmat) <- (IRanges::start(x = region)):(IRanges::end(x = region))
  
  ##########
  
  # find all the ATAC peak regions that intersect with region in any way
  all_atac_peak <- object[["ATAC"]]@ranges
  overlap_res <- GenomicRanges::findOverlaps(
    query = all_atac_peak,
    subject = region
  )
  region_gene_peaks <- all_atac_peak[overlap_res@from]
  for(i in 1:length(region_gene_peaks)){
    region_gene_peaks[i] <- intersect(x = region_gene_peaks[i],
                                      y = region)
  }
  
  # compute the cutmat for each of these regions
  len <- length(region_gene_peaks)
  cutmat_inpeak_list <- lapply(1:len, function(i){
    tmp <- Signac:::CutMatrix(
      object = object,
      region = region_gene_peaks[i],
      assay = assay,
      cells = cells,
      verbose = FALSE
    )
    colnames(tmp) <- (IRanges::start(x = region_gene_peaks[i])):(IRanges::end(x = region_gene_peaks[i]))
    tmp
  })
  
  ##########
  
  # compute sums for each cell
  cell_total <- Matrix::rowSums(cutmat)
  cell_total_inpeak <- sapply(1:len, function(i){
    Matrix::rowSums(cutmat_inpeak_list[[i]])
  })
  res <- cell_total - Matrix::rowSums(cell_total_inpeak)
  
  res
}
  