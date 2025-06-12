.assign_future_cells <- function(seurat_obj,
                                 future_ident, # the varible name inside ident
                                 ident, # the column name in metadata
                                 lineage_future_size,
                                 lineage_variable){
  stopifnot(ident %in% colnames(seurat_obj@meta.data),
            future_ident %in% seurat_obj@meta.data[,ident],
            length(names(lineage_future_size)) > 0,
            length(names(lineage_future_size)) == length(lineage_future_size))
  
  total_future_cells <- sum(lineage_future_size)
  current_future_cells <- length(which(seurat_obj@meta.data[,ident] == future_ident))
  stopifnot(total_future_cells <= current_future_cells)
  
  if(total_future_cells < current_future_cells){
    seurat_obj <- .throwing_out_cells(seurat_obj,
                                      future_ident,
                                      ident,
                                      lineage_future_size)
  }
  
  # now assign cells
  seurat_obj <- .assigning_cells_helper(seurat_obj,
                                        future_ident,
                                        ident,
                                        lineage_future_size,
                                        lineage_variable)
  
  seurat_obj
}

.throwing_out_cells <- function(seurat_obj,
                                future_ident,
                                ident,
                                lineage_future_size){
  total_future_cells <- sum(lineage_future_size)
  current_future_cells <- length(which(seurat_obj@meta.data[,ident] == future_ident))
  stopifnot(total_future_cells < current_future_cells)
  
  idx <- which(seurat_obj@meta.data[,ident] == future_ident)
  
  # remove cells
  num_remove <- current_future_cells - total_future_cells
  if(num_remove > 0) idx_remove <- sample(idx, num_remove, replace = FALSE)
  keep_vec <- rep(TRUE, length(Seurat::Cells(seurat_obj)))
  keep_vec[idx_remove] <- FALSE
  seurat_obj$keep <- keep_vec
  seurat_obj <- subset(seurat_obj, keep == TRUE)
  
  seurat_obj
}

.assigning_cells_helper <- function(seurat_obj,
                                    future_ident,
                                    ident,
                                    lineage_future_size,
                                    lineage_variable){
  idx <- which(seurat_obj@meta.data[,ident] == future_ident)
  
  lineage_vec <- unlist(lapply(1:length(lineage_future_size), function(i){
    rep(names(lineage_future_size)[i], lineage_future_size[i])
  }))
  
  stopifnot(length(lineage_vec) == length(idx))
  
  lineage_vec <- sample(lineage_vec, replace = FALSE)
  seurat_obj@meta.data[idx,lineage_variable] <- lineage_vec
  
  seurat_obj
}