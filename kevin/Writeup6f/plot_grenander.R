plot_grenander <- function(obj, 
                           xlim = range(obj$x),
                           ylim = range(obj$pdf), 
                           main = "",
                           ...){
  bin_cutoff <- seq(min(obj$x), max(obj$x), length.out = 1000)[-c(1,1000)]
  obj2 <- multiomeFate:::.add_cutoffs_to_grenander(obj = obj,
                                                   bin_cutoff = bin_cutoff)
  area_vec <- cumsum(diff(obj2$x)*obj2$pdf[-length(obj2$pdf)])
  idx <- which.min(abs(area_vec - 0.5))
  midpoint <- obj2$x[idx]
  
  plot(NA, xlim = xlim, ylim = ylim, main = paste0(main, "\nMidpoint: ", round(midpoint,4)), ...)
  
  x_seq <- seq(xlim[1], xlim[2], length.out = 11)
  y_seq <- seq(ylim[1], ylim[2], length.out = 11)
  for(x in x_seq) lines(rep(x,2), ylim, lwd = 0.5, lty = 2)
  for(y in y_seq) lines(xlim, rep(y,2), lwd = 0.5, lty = 2)
  
  x <- rep(obj$x, each = 2)
  y <- rep(obj$pdf, each = 2)
  y <- c(y[length(y)],y[-length(y)])
  polygon(x = x, y = y, col = "gray")
  
  lines(rep(midpoint,2), ylim, lwd = 2, lty = 3, col = 2)
  
  invisible()
}
