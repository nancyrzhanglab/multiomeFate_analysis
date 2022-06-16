barplot_function <- function(lineage_name,
                             lineage_mat,
                             membership_vec,
                             umap_mat,
                             col_indicator = rgb(1, 0.1, 0.1, 0.2),
                             col_palette_enrichment = paste0(grDevices::colorRampPalette(c('lightgrey', 'blue'))(100), "33"),
                             mode = c("indicator", "enrichment")[1]){
  dataset_vec <- c("day0", "day10_CIS", "week5_CIS",
                   NA, "day10_COCL2", "week5_COCL2",
                   NA, "day10_DABTRAM", "week5_DABTRAM")
  stopifnot(all(sort(unique(membership_vec)) == sort(dataset_vec[!is.na(dataset_vec)])),
            length(mode) == 1, mode %in% c("indicator", "enrichment"))
  
  lineage_mat2 <- Matrix::t(lineage_mat)
  j <- which(colnames(lineage_mat2) == lineage_name)
  cell_idx <- .nonzero_col(lineage_mat2, col_idx = j, bool_value = F)
  
  par(mfrow = c(3,3), mar = c(4,4,4,0.5))
  for(dataset in dataset_vec){
    cell_idx_all <- which(membership_vec == dataset)
    
    if(is.na(dataset)) {
      plot(NA, xaxt = "n", yaxt = "n", xlim = c(0,1), ylim = c(0,1), 
           bty = "n", main = "", xlab = "", ylab = "")
    } else {
      cell_idx2 <- intersect(cell_idx, cell_idx_all)
      
      plot(x = umap_mat[,1], y = umap_mat[,2],
           xlab = colnames(umap_mat)[1], ylab = colnames(umap_mat)[2],
           xaxt = "n", yaxt = "n", bty = "n",
           pch = 16, col = "gray",
           main = paste0(lineage_name, " for ", dataset, "\n(Present in ", 
                         length(cell_idx2), " of ", length(cell_idx_all), " cells)"))
      axis(1); axis(2)
      
      if(length(cell_idx2) > 0){
        if(mode == "indicator"){
          points(x = umap_mat[cell_idx2,1], y = umap_mat[cell_idx2,2],
                 pch = 16, col = "white", cex = 2)
          points(x = umap_mat[cell_idx2,1], y = umap_mat[cell_idx2,2],
                 pch = 16, col = col_indicator, cex = 1.5)
          
        } else if(mode == "enrichment"){
          cell_enrichment <- sapply(cell_idx2, function(i){
            tmp_idx <- .nonzero_col(lineage_mat, col_idx = i, bool_value = F)
            tmp_vec <- .nonzero_col(lineage_mat, col_idx = i, bool_value = T)
            stopifnot(j %in% tmp_idx)
            if(length(tmp_idx) == 1){
              denominator <- tmp_vec[1]
            } else {
              denominator <- sum(sort(tmp_vec, decreasing = T)[1:2])
            }
            
            tmp_val <- tmp_vec[which(tmp_idx == j)]/denominator
          })
          cell_col <- sapply(cell_enrichment, function(val){
            col_palette_enrichment[ceiling(val*length(col_palette_enrichment))]
          })
          
          points(x = umap_mat[cell_idx2,1], y = umap_mat[cell_idx2,2],
                 pch = 16, col = "white", cex = 2)
          points(x = umap_mat[cell_idx2,1], y = umap_mat[cell_idx2,2],
                 pch = 16, col = cell_col, cex = 1.5)
        }
      }
    }
  }
  
  invisible()
}

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}
