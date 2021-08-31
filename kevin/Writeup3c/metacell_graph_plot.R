plot_metacell_graph <- function(embedding, clustering, adj_mat,
                                feature_vec = NA,
                                zlim = range(feature_vec),
                                base_color = grDevices::rgb(0.803, 0.156, 0.211),
                                col_vec = scales::hue_pal()(length(unique(clustering))),
                                asp = T,
                                bins = 100,
                                ...){
  if(all(!is.na(feature_vec))){
    stopifnot(length(feature_vec) == nrow(adj_mat))
    col_ramp <- grDevices::colorRampPalette(c("white", base_color))(bins)
    val_vec <- seq(zlim[1], zlim[2], length = bins)
    col_idx <- sapply(feature_vec, function(x){
      which.min(abs(x - val_vec))
    })
    col_vec <- col_ramp[col_idx]
  }
  
  median_coord <- compute_median_coords(embedding, clustering)
  
  for(i in 1:3){
    graphics::plot(NA, 
                   xlim = range(embedding[,1]),
                   ylim = range(embedding[,2]),
                   asp = asp,
                   ...)
    
    for(j in 1:nrow(clust_mat)){
      idx <- which(adj_mat[j,] != 0)
      for(j2 in idx){
        lines(median_coord[c(j,j2),])
      }
    }
    
    for(k in 1:length(uniq_clust)){
      points(median_coord[k,1], median_coord[k,2], 
             pch = 16, 
             cex = 2, 
             col = col_vec[k])
    }
  }
  
  invisible()
}

# https://stackoverflow.com/questions/13355176/gradient-legend-in-base
plot_legend <- function(zlim,
                        base_color = grDevices::rgb(0.803, 0.156, 0.211),
                        bins = 100, 
                        legend_coords = c(0,0,1,1),
                        spacing = 5,
                        main = "Legend",
                        xlim = c(0,2), 
                        ylim = c(0,1),
                        offset = 0.5,
                        cex = 1,
                        ...){
  col_ramp <- grDevices::colorRampPalette(c("white", base_color))(bins)
  
  legend_image <- as.raster(matrix(col_ramp, ncol=1))
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       type = 'n', 
       axes = F,
       xlab = '',
       ylab = '', 
       main = main,
       ...)
  text(x = xlim[2]+offset, 
       y = seq(ylim[1], ylim[2], l=spacing), 
       labels = seq(zlim[1], zlim[2], l=spacing),
       cex = cex)
  rasterImage(legend_image, 
              legend_coords[1], 
              legend_coords[2], 
              legend_coords[3], 
              legend_coords[4])
  
  invisible()
}

##################################

compute_median_coords <- function(embedding, clustering){
  stopifnot(is.character(clustering), length(clustering) == nrow(embedding))
  uniq_clust <- sort(unique(clustering))
  
  median_coord <- t(sapply(uniq_clust, function(clust){
    idx <- which(clustering == clust)
    apply(embedding[idx,,drop = F], 2, median)
  }))
  rownames(median_coord) <- clustering
  
  median_coord
}