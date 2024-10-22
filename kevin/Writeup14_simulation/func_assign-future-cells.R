.assign_future_cells <- function(seurat_obj,
                                 future_ident,
                                 ident,
                                 lineage_future_size){
  stopifnot(ident %in% colnames(seurat_obj@meta.data),
            future_ident %in% seurat_obj@meta.data[,ident],
            length(names(lineage_future_size)) > 0,
            length(names(lineage_future_size)) == length(lineage_future_size))
  
  total_future_cells <- sum(lineage_future_size)
  current_future_cells <- length(which(seurat_obj@meta.data[,ident] == future_ident))
  
  if(total_future_cells < current_future_cells){
    seurat_obj <- .assign_future_cells_throwing_out(seurat_obj,
                                                    future_ident,
                                                    ident,
                                                    lineage_future_size)
  } else {
    
  }
  
  # now assign cells
  
  seurat_obj
}

.assign_future_cells_throwing_out <- function(seurat_obj,
                                              future_ident,
                                              ident,
                                              lineage_future_size){
  total_future_cells <- sum(lineage_future_size)
  current_future_cells <- length(which(seurat_obj@meta.data[,ident] == future_ident))
  stopifnot(total_future_cells < current_future_cells)
  
  idx <- which(seurat_obj@meta.data[,ident] == future_ident)
  
  # remove cells
  num_remove <- current_future_cells - total_future_cells
  idx_remove <- sample(idx, num_remove, replace = FALSE)
  keep_vec <- rep(TRUE, length(Seurat::Cells(seurat_obj)))
  keep_vec[idx_remove] <- FALSE
  seurat_obj$keep <- keep_vec
  seurat_obj <- subset(seurat_obj, keep == TRUE)
  
  seurat_obj
}