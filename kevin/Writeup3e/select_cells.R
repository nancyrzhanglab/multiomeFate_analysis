.tabulate_lineages <- function(seurat_obj, 
                               condition_var = "Original_condition"){
  lin_mat <- seurat_obj[["lineage"]]@counts
  tmp <- Matrix::t(lin_mat)
  lin_idx_list <- lapply(1:ncol(tmp), function(j){
    .nonzero_col(tmp, j)
  })
  names(lin_idx_list) <- colnames(tmp)
  factor_vec <- as.factor(seurat_obj@meta.data[,condition_var])
  tabulate_mat <- t(sapply(lin_idx_list, function(idx){
    table(factor_vec[idx])
  }))
  rownames(tabulate_mat) <- colnames(tmp)
  tabulate_mat <-  tabulate_mat[order(tabulate_mat[,"naive"], decreasing = T),]
  
  tabulate_mat
}

########################3

.select_cells_per_lineage <- function(seurat_obj,
                                      lineage_selected){
  lin_mat <- Matrix::t(seurat_obj[["lineage"]]@counts)
  all_cell_names <- unique(unlist(lapply(lineage_selected, function(lineage){
    j <- which(colnames(lin_mat) == lineage)
    rownames(lin_mat)[.nonzero_col(lin_mat, j)]
  })))
  
  all_cell_names
}

.select_expansion_naives <- function(seurat_obj,
                                     tabulate_mat,
                                     treatment, 
                                     threshold = 2,
                                     type = 1){
  stopifnot(treatment != "naive", treatment %in% colnames(tabulate_mat))
  
  tmp <- tabulate_mat[,treatment]/tabulate_mat[,"naive"]
  tmp[is.infinite(tmp)] <- 0
  lineage_idx <- which(tmp > threshold)
  length(lineage_idx)
  lineage_selected <- rownames(tabulate_mat)[lineage_idx]
  
  # [[note to self: ugly coding -- reformat this]]
  if(type == 1){
    # selecting cells
    celltypes <- unique(seurat_obj@meta.data$Original_condition)
    cell_terminal <- which(seurat_obj@meta.data$Original_condition == treatment)
    
    all_cell_names <- .select_cells_per_lineage(seurat_obj,
                                                lineage_selected)
    naive_all <- rownames(seurat_obj@meta.data)[which(seurat_obj@meta.data$Original_condition == "naive")]
    naive_terminal <- all_cell_names[which(all_cell_names %in% naive_all)]
    
    list(lineage_selected = lineage_selected,
         naive_all = naive_all,
         naive_terminal = naive_terminal)
  } else {
    # selected collapsed lineages
    naive_idx <- which(seurat_obj@meta.data[,"Original_condition"] == "naive")
    idx <- intersect(which(seurat_obj@meta.data[,"lineage"] %in% lineage_selected),
                     naive_idx)
    
    naive_all <- rownames(seurat_obj@meta.data)[naive_idx]
    naive_terminal <- rownames(seurat_obj@meta.data)[idx]
  }
  
  list(lineage_selected = lineage_selected,
       naive_all = naive_all,
       naive_terminal = naive_terminal)
}

.select_naive_deathtoall <- function(seurat_obj,
                                     tabulate_mat, 
                                     threshold = 2,
                                     type = 1){
  stopifnot(treatment != "naive")
  naive_idx <- which(colnames(tabulate_mat) == "naive")
  
  bool_vec <- sapply(1:nrow(tabulate_mat), function(i){
    all(tabulate_mat[i,-naive_idx] <= threshold * tabulate_mat[i,naive_idx])
  })
  
  lineage_selected <- rownames(tabulate_mat)[which(bool_vec)]
  
  if(type == 1){
    all_cell_names <- .select_cells_per_lineage(seurat_obj,
                                                lineage_selected)
  } else {
    all_cell_names <- rownames(seurat_obj@meta.data)[which(seurat_obj@meta.data[,"lineage"] %in% lineage_selected)]
  }
  
  naive_all <- rownames(seurat_obj@meta.data)[which(seurat_obj@meta.data$Original_condition == "naive")]
  
  all_cell_names[which(all_cell_names %in% naive_all)]
  
}

##############################3

.collapse_by_lineage <- function(seurat_obj, assay, slot){
  lin_mat <- Matrix::t(seurat_obj[["lineage"]]@counts)
  cell_names_list <- lapply(colnames(lin_mat), function(lineage){
    j <- which(colnames(lin_mat) == lineage)
    .nonzero_col(lin_mat, j)
  })
  names(cell_names_list) <- colnames(lin_mat)
  
  # should probably just add this into the seurat_obj
  lineage_vec <- rep(NA, nrow(seurat_obj@meta.data))
  for(i in 1:length(cell_names_list)){
    lineage_vec[cell_names_list[[i]]] <- names(cell_names_list)[i]
  }
  
  lineage_type_vec <- sapply(1:length(lineage_vec), function(i){
    paste0(lineage_vec[i], "-", seurat_obj@meta.data$Original_condition[i])
  })
  unique_pairings <- unique(lineage_type_vec)
  count_mat <- SeuratObject:::GetAssayData.Seurat(seurat_obj, slot = slot, assay = assay)
  
  if(class(count_mat) != "matrix") count_mat <- as.matrix(count_mat)
  print("here")
  new_count_mat <- sapply(unique_pairings, function(pairing){
    idx <- which(lineage_type_vec == pairing)
    Matrix::rowMeans(count_mat[,idx,drop = F])
  })
  rownames(new_count_mat) <- rownames(count_mat)
  colnames(new_count_mat) <- unique_pairings
  
  new_seurat_obj <- Seurat::CreateSeuratObject(counts = new_count_mat)
  new_seurat_obj[["lineage"]] <- sapply(strsplit(unique_pairings, split = "-"), function(x){x[1]})
  new_seurat_obj[["Original_condition"]] <- sapply(strsplit(unique_pairings, split = "-"), function(x){x[2]})
  new_seurat_obj[["num_cells"]] <- sapply(unique_pairings, function(pairing){length(which(lineage_type_vec == pairing))})
  
  new_seurat_obj
}

########################3

.order_genes_by_threshold <- function(genes, 
                                      naive_terminal,
                                      seurat_obj,
                                      assay){
  cell_idx <- which(seurat_obj@meta.data$Original_condition == "naive")
  stopifnot(all(genes %in% rownames(seurat_obj[[assay]]@scale.data)))
  mat <- seurat_obj[[assay]]@scale.data[genes, cell_idx]
  
  cells.1 <- which(colnames(mat) %in% naive_terminal)
  cells.2 <- setdiff(1:ncol(mat), cells.1)
  
  proportion_list <- lapply(1:nrow(mat), function(j){
    vec <- mat[j,]
    min_val <- max(vec[vec <= 0])
    max_val <- min(vec[vec >= 0])
    sd1 <- stats::sd(vec[vec <= 0])
    sd2 <- stats::sd(vec[vec > 0])
    sd_val <- min(c(sd1, sd2))
    
    bool <- abs(min_val - max_val) >= 3*sd_val
    if(!bool || sd1 > sd2) return(list(mat = NA, pval = 1))
    
    res_mat <- matrix(0, 2, 2)
    colnames(res_mat) <- c("survive", "die")
    rownames(res_mat) <- c("low", "high")
    res_mat[1,1] <- length(which(vec[cells.1] <= 0))
    res_mat[2,1] <- length(which(vec[cells.1] > 0))
    
    res_mat[1,2] <- length(which(vec[cells.2] <= 0))
    res_mat[2,2] <- length(which(vec[cells.2] > 0))
    
    prop_test <- stats::prop.test(x = res_mat[2,], n = colSums(res_mat))
  
    list(mat = res_mat, pval = prop_test$p.value)
  })
  
  order_idx <- order(sapply(proportion_list, function(x){
    if(!all(is.na(x$mat))){
      max(x$mat[2,1]/x$mat[2,2], x$mat[2,2]/x$mat[2,1])
    } else {
      0
    }
  }), decreasing = T)
  names(proportion_list) <- rownames(mat)
  proportion_list <- proportion_list[order_idx]

  
  proportion_list
}