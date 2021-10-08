rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/kevin/Writeup3d/09302021_sydney_de.RData")
ls()

# find all the relevant terminal cells in preparation for a differential expression
celltypes <- unique(all_data@meta.data$Original_condition)
terminal_celltypes <- setdiff(celltypes, "naive")
de_terminal_list <- lapply(terminal_celltypes, function(celltype){
  print(celltype)
  cells.1 <- which(all_data@meta.data$Original_condition == celltype)
  cells.2 <- which(all_data@meta.data$Original_condition %in% setdiff(terminal_celltypes, celltype))
  
  set.seed(10)
  Seurat:::FindMarkers(all_data[["SCT"]],
                       slot = "data",
                       cells.1 = rownames(all_data@meta.data)[cells.1],
                       cells.2 = rownames(all_data@meta.data)[cells.2],
                       verbose = T)
})

for(i in de_terminal_list){
  print(head(i))
  print(tail(i))
  print(quantile(i$p_val_adj))
  print("===")
}

#########################

de_naive_list <- lapply(terminal_celltypes, function(celltype){
  print(celltype)
  
  cells <- which(all_data@meta.data$Original_condition == celltype)
  lineage_idx <- which(sparseMatrixStats::rowSums2(all_data[["lineage"]]@counts[,cells]) > 0)
  lineages <- rownames(all_data[["lineage"]]@counts)[lineage_idx]
  
  lineage_cell_idx1 <- which(sparseMatrixStats::colSums2(all_data[["lineage"]]@counts[lineages,]) > 0)
  lineage_cell_idx2 <- which(all_data@meta.data$Original_condition == "naive")
  cells.1 <- intersect(lineage_cell_idx1, lineage_cell_idx2)
  cells.2 <- setdiff(lineage_cell_idx2, lineage_cell_idx1)
  print(length(cells.1))
  print(length(cells.2))
  
  set.seed(10)
  Seurat:::FindMarkers(all_data[["SCT"]],
                       slot = "data",
                       cells.1 = rownames(all_data@meta.data)[cells.1],
                       cells.2 = rownames(all_data@meta.data)[cells.2],
                       logfc.threshold = 0.1,
                       verbose = T)
})

for(i in de_naive_list){
  print(dim(i))
  print(head(i))
  print(tail(i))
  print(quantile(i$p_val_adj))
  print("===")
}

save(de_terminal_list , de_naive_list,
     file = "../../../../out/kevin/Writeup3d/09302021_sydney_de_list.RData")

##############################

terminal_gene_list <- lapply(de_terminal_list, function(x){
  rownames(x)[1:50]
})
tab <- table(unlist(terminal_gene_list))
rm_genes <- names(tab)[which(tab > 1)]
terminal_gene_list <- lapply(terminal_gene_list, function(vec){
  vec <- setdiff(vec, rm_genes)
})

naive_gene_list <- lapply(de_naive_list, function(x){
  rownames(x)[1:50]
})
tab <- table(unlist(naive_gene_list))
rm_genes <- names(tab)[which(tab > 1)]
naive_gene_list <- lapply(naive_gene_list, function(vec){
  vec <- setdiff(vec, rm_genes)
})

length(intersect(unlist(naive_gene_list), unlist(terminal_gene_list)))
rm_genes <- intersect(unlist(naive_gene_list), unlist(terminal_gene_list))

all_gene_list <- vector("list", 2*length(terminal_celltypes))
for(i in 1:length(terminal_celltypes)){
  tmp <- naive_gene_list[[i]]
  tmp <- setdiff(tmp, rm_genes)
  all_gene_list[[2*i+1]] <- tmp
  
  tmp <- terminal_gene_list[[i]]
  tmp <- setdiff(tmp, rm_genes)
  all_gene_list[[2*i+2]] <- tmp
}
table(table(unlist(all_gene_list)))

all_cell_list <- vector("list", 2*length(terminal_celltypes))
for(i in 1:length(terminal_celltypes)){
  
  cells <- which(all_data@meta.data$Original_condition == terminal_celltypes[i])
  lineage_idx <- which(sparseMatrixStats::rowSums2(all_data[["lineage"]]@counts[,cells]) > 0)
  lineages <- rownames(all_data[["lineage"]]@counts)[lineage_idx]
  
  lineage_cell_idx1 <- which(sparseMatrixStats::colSums2(all_data[["lineage"]]@counts[lineages,]) > 0)
  lineage_cell_idx2 <- which(all_data@meta.data$Original_condition == "naive")
  cells.1 <- intersect(lineage_cell_idx1, lineage_cell_idx2)
  
  all_cell_list[[2*i+1]] <- cells.1
  all_cell_list[[2*i+2]] <- cells
}
table(table(unlist(all_cell_list)))

mat <- all_data[["SCT"]]@data[unlist(all_gene_list), unlist(all_cell_list)]
mat <- as.matrix(Matrix::t(mat))

mat <- scale(mat, center = T, scale = T)
graphics::image(.rotate(mat))

#################################

scaling_power <- 1
n <- nrow(mat)
zlim <- range(mat)
mat <- scale(mat, center = T, scale = T)

# construct colors. green is negative
max_val <- max(abs(zlim)); min_val <- max(min(abs(zlim)), 1e-3)
col_vec_neg <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(1,1,1), num_col,
                                 luminosity = luminosity)
break_vec_neg <- -max_val*seq(1, 0, length.out = num_col+2)^scaling_power
break_vec_neg <- break_vec_neg[-length(break_vec_neg)]

# red is positive
col_vec_pos <- .colorRamp_custom(c(1,1,1), c(0.803, 0.156, 0.211), num_col,
                                 luminosity = luminosity)
break_vec_pos <- max_val*seq(0, 1, length.out = num_col+2)^scaling_power
break_vec_pos <- break_vec_pos[-1]

# combine the two
break_vec <- c(break_vec_neg, break_vec_pos)
col_vec <- c(col_vec_neg, "white", col_vec_pos)

membership_vec <- as.numeric(membership_vec) ## convert to integers
idx <- order(membership_vec, decreasing = F)
breakpoints <- 1-which(abs(diff(sort(membership_vec, decreasing = F))) >= 1e-6)/n

line_func <- function(){
  if(!all(is.na(membership_vec))){
    for(i in 1:length(breakpoints)){
      graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2.1, col = "white")
      if(!all(is.na(major_breakpoint)) && i %in% major_breakpoint){
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2.2, lty = 1)
      } else {
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2, lty = 2)
      }
    }
  }
}

graphics::image(.rotate(score_mat[idx,,drop = F]), col = col_vec,
                breaks = break_vec, ...)
line_func()

invisible()

