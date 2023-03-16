plot_coveragetracks <- function(
    cutmat,
    peaks,
    cell_nopeak_offset = 1/100,
    col_nopeak_cell = 2,
    col_peak = rgb(255, 102, 102, alpha = 0.2*255, maxColorValue = 255),
    col_track = rgb(0.5,0.5,0.5,0.8),
    max_height = max(cutmat@x),
    main = "",
    tol = 1e-5
){
  n <- nrow(cutmat)
  p <- ncol(cutmat)
  x_vec <- as.numeric(colnames(cutmat))
  
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
  
  # mark cells that had no peak
  sum_vec <- Matrix::rowSums(cutmat)
  if(any(sum_vec <= tol)){
    cell_idx <- which(sum_vec <= tol)
    x_val <- x_vec[1] + cell_nopeak_offset*diff(range(range(x_vec)))
    .draw_circle(
      x_vec = rep(x_val, length(cell_idx)),
      y_vec = cell_idx-1,
      radius = 0.4,
      col = col_nopeak_cell
    )
  }
  
  invisible()
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