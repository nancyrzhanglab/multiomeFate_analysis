plot_coveragetracks <- function(
    cutmat,
    peaks,
    atac_count = NULL,
    atac_count_max = NULL,
    cell_nopeak_offset = 1/100,
    col_nopeak_cell = 2,
    col_peak = rgb(255, 102, 102, alpha = 0.2*255, maxColorValue = 255),
    col_scale = c('lightgrey', 'blue'),
    col_track = rgb(0.5,0.5,0.5,0.8),
    max_height = max(cutmat@x),
    main = "",
    rna_gene_count = NULL,
    rna_gene_count_max = NULL,
    tol = 1e-5
){
  n <- nrow(cutmat)
  p <- ncol(cutmat)
  x_vec <- as.numeric(colnames(cutmat))
  
  sum_vec <- Matrix::rowSums(cutmat)
  main <- paste0(main, " (", round(length(which(sum_vec > tol))/length(sum_vec)*100), "% non-zero chrom. frag)")
  graphics::plot(NA, 
                 xlab = "Basepair",
                 xlim = range(x_vec), 
                 ylab = "Coverage",
                 ylim = c(0,n),
                 main = main)
  
  # plot one line per cell
  for(i in 0:(n-1)){
    if(i %% 10 == 0) {lty <- 3} else if(i %% 5 == 0) {lty <- 2} else lty <- 1
    graphics::lines(
      x = range(x_vec),
      y = rep(i,2),
      lty = lty
    )
  }
  
  # plot coverages
  for(i in 1:n){
    graphics::polygon(
      x = c(x_vec[1], x_vec, x_vec[c(p,1)]),
      y = c(0,as.numeric(cutmat[i,]),0,0)/max_height+(i-1),
      border = NA,
      density = NULL,
      col = col_track
    )
  }
  
  # plot peak regions
  for(j in 1:nrow(peaks)){
    graphics::polygon(
      x = c(peaks[j,c(1,2,2,1)]),
      y = n*c(2,2,-2,-2),
      border = NA,
      density = NULL,
      col = col_peak
    )
  }
  
  xmult <- .getXmult()
  
  # mark cells that had no peak
  if(any(sum_vec <= tol)){
    cell_idx <- which(sum_vec <= tol)
    x_val <- x_vec[1] + 0.5*xmult
    .draw_circle(
      x_vec = rep(x_val, length(cell_idx)),
      y_vec = cell_idx-1,
      radius = 0.4,
      col = col_nopeak_cell
    )
  }
  
  # plot RNA count for this specific gene
  if(!any(is.null(rna_gene_count)) & !is.null(rna_gene_count_max)){
    stopifnot(length(rna_gene_count) == nrow(cutmat))
    rna_gene_count <- pmin(rna_gene_count, rna_gene_count_max)
    
    # compute the appropriate values
    col_pal <- grDevices::colorRampPalette(col_scale)(100)
    break_vec <- seq(0, rna_gene_count_max, length.out = 100)
    col_vec <- sapply(rna_gene_count, function(x){
      col_pal[which.min(abs(x - break_vec))]
    })
    height_pal <- seq(0, 1, length.out = 100)
    height_vec <- sapply(rna_gene_count, function(x){
      height_pal[which.min(abs(x - break_vec))]
    })
    
    # plot
    x_val <- x_vec[1] + 1*xmult
    for(i in 1:length(col_vec)){
      graphics::polygon(
        x = c(rep(x_val, 2), rep(x_val+1.5*xmult, 2)),
        y = (i-1)+c(0,rep(height_vec[i],2),0),
        border = NA,
        density = NULL,
        col = col_vec[i]
      )
    }
    
    graphics::text(x = x_val, y = n, labels = "RNA count", cex = 0.5)
  }
  
  # plot ATAC count
  if(!any(is.null(atac_count)) & !is.null(atac_count_max)){
    stopifnot(length(atac_count) == nrow(cutmat))
    atac_count <- pmin(atac_count, atac_count_max)
    
    # compute the appropriate values
    col_pal <- grDevices::colorRampPalette(col_scale)(100)
    break_vec <- seq(0, atac_count_max, length.out = 100)
    col_vec <- sapply(atac_count, function(x){
      col_pal[which.min(abs(x - break_vec))]
    })
    height_pal <- seq(0, 1, length.out = 100)
    height_vec <- sapply(atac_count, function(x){
      height_pal[which.min(abs(x - break_vec))]
    })
    
    # plot
    x_val <- max(x_vec) - xmult
    for(i in 1:length(col_vec)){
      graphics::polygon(
        x = c(rep(x_val, 2), rep(x_val+xmult, 2)),
        y = (i-1)+c(0,rep(height_vec[i],2),0),
        border = NA,
        density = NULL,
        col = col_vec[i]
      )
    }
    
    graphics::text(x = x_val, y = n, labels = "Total ATAC", cex = 0.5)
  }
  
  invisible()
}

rna_gene_count_extractor <- function(
    object,
    gene,
    assay = "RNA",
    cells = NULL,
    slot = "data",
    scale_factor = 1e6
){
  cells <- Signac:::SetIfNull(x = cells, y = colnames(x = object))
  cells <- intersect(cells, colnames(x = object))
  
  mat <- Seurat::GetAssayData(object = object, assay = assay, slot = slot)
  idx <- which(rownames(mat) == gene)
  if(length(idx) != 1) return(NULL)
  
  res <- as.numeric(mat[idx,cells])/object$nCount_RNA[cells]*scale_factor
  names(res) <- cells
  res
}

atac_count_extractor <- function(
    object,
    cells = NULL
){
  res <- object$nCount_ATAC
  
  cells <- Signac:::SetIfNull(x = cells, y = colnames(x = object))
  cells <- intersect(cells, colnames(x = object))
  
  res[cells]
}

######################################


# from https://github.com/plotrix/plotrix/blob/master/R/getYmult.R

.getXmult <- function() {
  if(dev.cur() == 1) {
    warning("No graphics device open.")
    ymult<-1
  }
  else {
    # get the plot aspect ratio
    xyasp <- par("pin") # width and height of the current plotting region
    # get the plot coordinate ratio
    xycr <- diff(par("usr"))[c(1,3)] # "width" and "height" of the plotting coordinates
    xmult <- xyasp[2]/xyasp[1]*xycr[1]/xycr[2]
  }
  return(xmult)
}

# from https://github.com/plotrix/plotrix/blob/master/R/draw.circle.R
.draw_circle <- function(x_vec,
                         y_vec,
                         radius, #radius is techincally radius in the Y direction
                         border = NA,
                         nv = 100,
                         col = "gray", 
                         lty = 1,
                         density = NULL,
                         angle = 45,
                         lwd = 1) {
  xylim <- par("usr")
  plotdim <- par("pin")
  xmult <- .getXmult()
  angle.inc <- 2*pi/nv
  angles <- seq(0,2*pi-angle.inc,by=angle.inc)
  
  if(length(col) < length(x_vec)) 
    col <- rep(col,length.out=length(x_vec))
  
  for(circle in 1:length(x_vec)) {
    xv <- cos(angles)*radius*xmult + x_vec[circle] # xmult is needed to make the circle look like a circle
    yv <- sin(angles)*radius + y_vec[circle]
    polygon(xv,yv,border=border,col=col[circle],lty=lty,
            density=density,angle=angle,lwd=lwd)
  }
  invisible()
}