plot_topic_heatmap <- function(topic_mat,
                               gene_vec,
                               file_name,
                               height_max = 6500,
                               width_base = 2000){
  # Custom heatmap
  col_vec <- viridis::viridis(25)
  break_vec <- seq(0, 1, length.out = 26)
  idx <- which(rownames(topic_mat) %in% gene_vec)
  topic_mat <- topic_mat[idx,]
  rowname_vec <- rownames(topic_mat)
  l1_vec <- apply(topic_mat, 1, sum)
  topic_mat <- diag(1/l1_vec) %*% topic_mat
  rownames(topic_mat) <- rowname_vec
  
  hclust_res <- stats::hclust(stats::dist(topic_mat))
  topic_mat <- topic_mat[hclust_res$order,]
  
  ratio <- min(width_base*(nrow(topic_mat)-1)/(ncol(topic_mat)-1), height_max)/width_base
  png(file_name,
      height = width_base*ratio, 
      width = width_base, 
      res = 500, 
      units = "px")
  
  y_indent_val <- 1/(2*(nrow(topic_mat)-1))
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  image(tiltedCCA:::.rotate(topic_mat), 
        asp = ratio, ylim = c(-y_indent_val,1+y_indent_val),
        main = "", xlab = "", ylab = "",
        xaxt = "n", yaxt = "n", bty = "n",
        breaks = break_vec, col = col_vec)
  
  # draw lines
  x_indent_val <- 1/(2*(ncol(topic_mat)-1))
  x_vec <- seq(-x_indent_val, 1+x_indent_val, by = 2*x_indent_val)
  for(x in x_vec[seq(ncol(topic_mat)+1,6,by=-5)[-1]]){
    graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 2, col = "white")
    graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 1.5, lty = 2)
  }
  
  y_vec <- seq(-y_indent_val, 1+y_indent_val, by = 2*y_indent_val)
  for(y in y_vec[seq(6,length(y_vec), by = 5)]){
    graphics::lines(x = c(0,1), y = rep(y, 2), lwd = 2, col = "white")
  }
  
  # label genes
  x_vec_label <- rep(x_vec[c(3,5)], times = ceiling(nrow(topic_mat)/2))[1:nrow(topic_mat)]
  y_vec_label <- seq(1, 0, by = -2*y_indent_val)
  graphics::text(x = x_vec_label, y = y_vec_label, labels = rownames(topic_mat),
                 col = "white", cex = 0.5)
  
  graphics.off()
}