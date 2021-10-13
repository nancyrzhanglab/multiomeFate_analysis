rm(list=ls())
library(Seurat); library(Signac)
source("funcs.R")

load("../../../../data/Sydney_stressors_2021-09-24/all_data_SCT.RData")

################3
all_data_safe <- all_data

# first clean up all the lineages
lin_mat <- all_data[["lineage"]]@counts
lin_mat@x <- rep(1, length(lin_mat@x))
dim(lin_mat)
lin_mat[1:5,1:5]
count_vec <- sparseMatrixStats::rowSums2(lin_mat)
lin_mat <- lin_mat[count_vec > 1, ] # since we don't have the barcodes themselves, the best we can do right now is focus on the lineages with many cells
sum_vec <- sparseMatrixStats::colSums2(lin_mat)
lin_mat <- lin_mat[,sum_vec > 0] 
dim(lin_mat)
keep_vec <- rep(0, ncol(lin_mat))
keep_vec[sum_vec > 0] <- 1
all_data[["keep"]] <- keep_vec
all_data <- subset(all_data, keep == 1)
all(colnames(all_data) == colnames(lin_mat))

# find the percentage of cells in each lineage
lin_mat <- Matrix::t(lin_mat)
lin_idx_list <- lapply(1:ncol(lin_mat), function(j){
  .nonzero_col(lin_mat, j)
})
names(lin_idx_list) <- colnames(lin_mat)
factor_vec <- as.factor(all_data@meta.data$Original_condition)
tabulate_mat <- t(sapply(lin_idx_list, function(idx){
  table(factor_vec[idx])
}))
rownames(tabulate_mat) <- colnames(lin_mat)
naive_idx <- which(colnames(tabulate_mat) == "naive")
tabulate_mat <- tabulate_mat[-which(tabulate_mat[,naive_idx] == 0),]
percentage_mat <- tabulate_mat
percentage_mat <- percentage_mat[rowSums(percentage_mat) != 0,]
percentage_mat <- .mult_vec_mat(1/rowSums(percentage_mat), percentage_mat)

# cluster the percentages
set.seed(10)
eps <- .l2norm(rep(1/6,6) - c(rep(1/5,5),0))/40
dbscan_clust <- dbscan::dbscan(percentage_mat, 
                               eps = eps,
                               minPts = 20)
table(dbscan_clust$cluster)

# find the mean of the percentages
cluster_vec <- factor(dbscan_clust$cluster)
names(cluster_vec) <- rownames(percentage_mat)
table(cluster_vec)
K <- length(levels(dbscan_clust$cluster))
percentage_mean <- t(sapply(levels(cluster_vec), function(k){
  idx <- which(cluster_vec == k)
  colMeans(percentage_mat[idx,,drop = F])
}))
round(percentage_mean, 2)

# go through the lineages -- assign each cell to a lineage
lin_mat <- lin_mat[,colnames(lin_mat) %in% rownames(percentage_mat)]
sum_vec <- sparseMatrixStats::rowSums2(lin_mat)
lin_mat <- lin_mat[sum_vec > 0,] 
lin_mat <- Matrix::t(lin_mat)
cell_idx_list <- lapply(1:ncol(lin_mat), function(j){
  .nonzero_col(lin_mat, j)
})
names(cell_idx_list) <- colnames(lin_mat)
stopifnot(all(sapply(cell_idx_list, length) > 0))

cell_idx_list2 <- cell_idx_list
lineage_coverage <- sparseMatrixStats::rowSums2(lin_mat) 
for(i in 1:length(cell_idx_list2)){
  if(i %% floor(length(cell_idx_list2)/10) == 0) cat('*')
  set.seed(10)
  vec <- as.numeric(as.character(cluster_vec))[cell_idx_list2[[i]]]
  vec <- factor(vec, levels = levels(cluster_vec))
  percentage_vec <- table(vec)/table(cluster_vec)
  
  max_cluster <- names(percentage_vec)[which.max(percentage_vec)]
  considered_lineages <- which(vec == max_cluster)
  max_idx <- which.max(lineage_coverage[cell_idx_list2[[i]][considered_lineages]])
  cell_idx_list2[[i]] <- (cell_idx_list2[[i]][considered_lineages])[max_idx]
}

#####################

i_vec <- unlist(cell_idx_list2)
j_vec <- 1:length(cell_idx_list2)
lin_mat2 <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = rep(1,length(i_vec)),
                                 dims = dim(lin_mat))
rownames(lin_mat2) <- rownames(lin_mat)
colnames(lin_mat2) <- colnames(lin_mat)

while(TRUE){
  print(dim(lin_mat2))
  
  bool <- TRUE
  sum_vec <- sparseMatrixStats::rowSums2(lin_mat2)
  if(any(sum_vec < 2)) bool <- FALSE
  lin_mat2 <- lin_mat2[sum_vec >= 2,] 
  
  sum_vec <- sparseMatrixStats::colSums2(lin_mat2)
  if(any(sum_vec == 0)) bool <- FALSE
  lin_mat2 <- lin_mat2[,sum_vec > 0] 
  
  if(bool) break()
}

cluster_vec2 <- cluster_vec[which(names(cluster_vec) %in% rownames(lin_mat2))]
table(cluster_vec2)
quantile(sparseMatrixStats::rowSums2(lin_mat2))

############################

lin_mat2 <- Matrix::t(lin_mat2)
lin_idx_list <- lapply(1:ncol(lin_mat2), function(j){
  .nonzero_col(lin_mat2, j)
})
names(lin_idx_list) <- colnames(lin_mat2)
factor_vec <- as.factor(all_data@meta.data$Original_condition)[rownames(all_data@meta.data) %in% rownames(lin_mat2)]
tabulate_mat <- t(sapply(lin_idx_list, function(idx){
  table(factor_vec[idx])
}))
rownames(tabulate_mat) <- colnames(lin_mat2)
naive_idx <- which(colnames(tabulate_mat) == "naive")
tabulate_mat <- tabulate_mat[-which(tabulate_mat[,naive_idx] == 0),]
lin_mat2 <- lin_mat2[,colnames(lin_mat2) %in% rownames(tabulate_mat)]

lin_mat2 <- Matrix::t(lin_mat2)
while(TRUE){
  print(dim(lin_mat2))
  
  bool <- TRUE
  sum_vec <- sparseMatrixStats::rowSums2(lin_mat2)
  if(any(sum_vec < 2)) bool <- FALSE
  lin_mat2 <- lin_mat2[sum_vec >= 2,] 
  
  sum_vec <- sparseMatrixStats::colSums2(lin_mat2)
  if(any(sum_vec == 0)) bool <- FALSE
  lin_mat2 <- lin_mat2[,sum_vec > 0] 
  
  if(bool) break()
}

cluster_vec2 <- cluster_vec[which(names(cluster_vec) %in% rownames(lin_mat2))]
table(cluster_vec2)


#############################

# now for all the checks:
# check 1: each lineage has at least 2 cells
sum_vec <- sparseMatrixStats::rowSums2(lin_mat2)
stopifnot(all(sum_vec >= 2))

# check 2: each cell is only in one lineage
sum_vec <- sparseMatrixStats::colSums2(lin_mat2)
stopifnot(all(sum_vec == 1))

# check 3: each lineage has at least one naive cell
tmp <- lin_mat2
tmp <- Matrix::t(tmp)
lin_idx_list <- lapply(1:ncol(tmp), function(j){
  .nonzero_col(tmp, j)
})
bool_vec <- sapply(lin_idx_list, function(idx){
  cell_names <- rownames(tmp)[idx]
  type_vec <- all_data@meta.data$Original_condition[which(rownames(all_data@meta.data) %in% cell_names)]
  any(type_vec == "naive")
})
stopifnot(all(bool_vec))


keep_vec <- rep(0, ncol(all_data[["lineage"]]@counts))
keep_vec[which(rownames(all_data@meta.data) %in% colnames(lin_mat2))] <- 1
all_data[["keep"]] <- keep_vec
all_data <- subset(all_data, keep == 1)
all_data[["lineage"]]@counts <- lin_mat2
all_data[["lineage"]]@data <- lin_mat2

table(all_data@meta.data$Original_condition)

############################

save(all_data,
     file = "../../../../out/kevin/Writeup3d/09302021_sydney_de_preprocess.RData")
