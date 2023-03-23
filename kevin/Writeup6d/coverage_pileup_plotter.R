coverage_pileup_plotter <- function(
    pileup_res, # output from compute_pileup_curve(),
    col_curve = 1,
    col_peak_median = rgb(255, 102, 102, alpha = 0.3*255, maxColorValue = 255),
    col_peak_max = rgb(240, 230, 140, alpha = 0.3*255, maxColorValue = 255),
    main = "",
    ylim = NA,
    ... # additional parameters for the curve
){
  x_vec <- as.numeric(names(pileup_res$pileup_vec))
  y_vec <- pileup_res$pileup_vec
  inflation_factor <- pileup_res$inflation_factor
  peak_width_max <- pileup_res$peak_width_max
  peak_width_median <- pileup_res$peak_width_median
  
  xlim <- range(x_vec)
  if(all(is.na(ylim))) ylim <- range(y_vec)
  
  plot(NA, 
       xlim = xlim, 
       ylim = ylim,
       xlab = "Basepair dist from peak",
       ylab = "Normalized value",
       main = paste0(main, "\nScore: ",
                     round(pileup_res$score,2)))
  
  # plot peak max region
  if(!is.null(peak_width_max)){
    window <- c(-1,1)*inflation_factor*peak_width_max
    graphics::polygon(x = c(rep(window[1],2), rep(window[2], 2)),
                      y = c(-1e5,2*ylim[2],2*ylim[2],-1e5),
                      border = NA,
                      density = NULL,
                      col = col_peak_max)
  }
  
  # plot peak median region
  if(!is.null(peak_width_median)){
    window <- c(-.5,.5)*peak_width_median
    graphics::polygon(x = c(rep(window[1],2), rep(window[2], 2)),
                      y = c(-1e5,2*ylim[2],2*ylim[2],-1e5),
                      border = NA,
                      density = NULL,
                      col = col_peak_median)
  }
  
  # plot curve
  lines(x = x_vec,
        y = y_vec,
        col = col_curve,
        ...)
  
  invisible()
}