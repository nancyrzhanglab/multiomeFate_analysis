.normalize_fragments <- function(cutmat,
                                 peak_mat,
                                 normalize_by_unique_cells = T){
  stopifnot(length(rownames(cutmat)) > 0)
  
  peak_width <- stats::median(apply(peak_mat, 1, diff))
  peak_loc <- round(apply(peak_mat, 1, mean))
  
  frag_locations <- multiomeFate:::.extract_fragment_from_cutmat(cutmat)
  cell_names <- multiomeFate:::.extract_fragment_cell_from_cutmat(cutmat)
  data_mat <- data.frame(
    cell_name = cell_names,
    frag_location = frag_locations
  )
  tab_vec <- table(data_mat$cell_name)
  
  # tabulate how many unique cells have a fragment nearby
  n <- nrow(data_mat)
  nearby_vec <- sapply(1:n, function(i){
    loc <- data_mat[i,"frag_location"]
    start <- loc - peak_width
    end <- loc + peak_width
    idx <- intersect(
      which(data_mat[,"frag_location"] >= start),
      which(data_mat[,"frag_location"] <= end)
    )
    
    ifelse(normalize_by_unique_cells, length(unique(data_mat$cell_name[idx])), length(idx)) 
  })
  
  data_mat$nearby_frags <- nearby_vec
  
  # assign fragments to peaks and compute the distance
  assigned_peak <- sapply(1:n, function(i){
    idx <- which.min(abs(peak_loc - data_mat$frag_location[i]))
    peak_loc[idx]
  })
  data_mat$assigned_peak <- assigned_peak
  data_mat$dist_to_peak <- pmax(abs(data_mat$frag_location - data_mat$assigned_peak), 1)
  
  # now do a tf-idf normalization
  tf <- sapply(1:n, function(i){
    1/tab_vec[data_mat$cell_name[i]]
  })
  if(normalize_by_unique_cells){
    idf <- log(nrow(cutmat)/data_mat$nearby_frags)
  } else {
    idf <- log(nrow(data_mat)/data_mat$nearby_frags)
  }
  
  normalized_value <- tf*idf
  data_mat$tfidf <- normalized_value
  
  data_mat
}


