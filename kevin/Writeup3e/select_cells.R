.tabulate_lineages <- function(seurat_obj){
  lin_mat <- seurat_obj[["lineage"]]@counts
  tmp <- Matrix::t(lin_mat)
  lin_idx_list <- lapply(1:ncol(tmp), function(j){
    .nonzero_col(tmp, j)
  })
  names(lin_idx_list) <- colnames(tmp)
  factor_vec <- as.factor(seurat_obj@meta.data$Original_condition)
  tabulate_mat <- t(sapply(lin_idx_list, function(idx){
    table(factor_vec[idx])
  }))
  rownames(tabulate_mat) <- colnames(tmp)
  tabulate_mat <-  tabulate_mat[order(tabulate_mat[,"naive"], decreasing = T),]
  
  tabulate_mat
}

.select_expansion_naives <- function(seurat_obj,
                                     tabulate_mat,
                                     treatment, 
                                     threshold = 2){
  stopifnot(treatment != "naive", treatment %in% colnames(tabulate_mat))
  
  tmp <- tabulate_mat[,treatment]/tabulate_mat[,"naive"]
  tmp[is.infinite(tmp)] <- 0
  lineage_idx <- which(tmp > threshold)
  length(lineage_idx)
  lineage_selected <-  rownames(tabulate_mat)[lineage_idx]
  
  celltypes <- unique(seurat_obj@meta.data$Original_condition)
  cell_terminal <- which(seurat_obj@meta.data$Original_condition == treatment)
  
  lin_mat <- Matrix::t(seurat_obj[["lineage"]]@counts)
  all_cell_names <- unique(unlist(lapply(lineage_selected, function(lineage){
    j <- which(colnames(lin_mat) == lineage)
    rownames(lin_mat)[.nonzero_col(lin_mat, j)]
  })))
  naive_all <- rownames(seurat_obj@meta.data)[which(seurat_obj@meta.data$Original_condition == "naive")]
  
  naive_terminal <- all_cell_names[which(all_cell_names %in% naive_all)]
  
  list(lineage_selected = lineage_selected,
       naive_all = naive_all,
       naive_terminal = naive_terminal)
}

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