rm(list=ls())
library(Seurat); library(Signac)
source("funcs.R")

load("../../../../data/Sydney_stressors_2021-09-24/all_data_SCT.RData")

################3

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
idx_list <- lapply(1:ncol(lin_mat), function(j){
  .nonzero_col(lin_mat, j)
})
names(idx_list) <- colnames(lin_mat)
factor_vec <- as.factor(all_data@meta.data$Original_condition)
tabulate_mat <- t(sapply(idx_list, function(idx){
  table(factor_vec[idx])
}))
rownames(tabulate_mat) <- colnames(lin_mat)
naive_idx <- which(colnames(tabulate_mat) == "naive")
tabulate_mat <- tabulate_mat[-which(tabulate_mat[,naive_idx] == 0),]
percentage_mat <- tabulate_mat[,-naive_idx]
percentage_mat <- percentage_mat[rowSums(percentage_mat) != 0,]
percentage_mat <- .mult_vec_mat(1/rowSums(percentage_mat), percentage_mat)
for(i in 1:ncol(percentage_mat)){
  print(colnames(percentage_mat)[i])
  print(paste0("Percentage of non-zeros: ", round(sum(percentage_mat[,i] != 0)/nrow(percentage_mat), 2)))
  print(paste0("Quantiles: "))
  print(quantile(percentage_mat[which(percentage_mat[,i] != 0),i]))
  print(paste0("Quantiles: "))
  print(quantile(percentage_mat[which(percentage_mat[,i] != 0),i]))
  print("====")
}

# keep only lineages that have concentrate mainly in one terminal-condition
keep_lin_list <- lapply(1:ncol(percentage_mat), function(k){
  rownames(percentage_mat)[which(percentage_mat[,k] > 0.9)]
})
names(keep_lin_list) <- colnames(percentage_mat)
sapply(keep_lin_list, length)

# now find the cells associated with these lineages
all_lin <- unlist(keep_lin_list)
lin_mat <- lin_mat[,all_lin]
dim(lin_mat)
# clean up lin_mat by zeroing out all the unnecessary cells
for(k in 1:length(keep_lin_list)){
  for(lin in keep_lin_list[[k]]){
    j <- which(colnames(lin_mat) == lin)
    cell_idx <- .nonzero_col(lin_mat, j)
    bad_idx <- which(!all_data@meta.data$Original_condition[cell_idx] %in% c("naive", names(keep_lin_list)[k]))
    if(length(bad_idx) > 0){
      lin_mat[cell_idx[bad_idx], j] <- 0
    }
  }
}

keep_vec <- rep(0, nrow(lin_mat))
keep_vec[sparseMatrixStats::rowSums2(lin_mat) > 0] <- 1
table(keep_vec)

all_data[["lineage"]]@data <- Matrix::t(lin_mat)
all_data[["lineage"]]@counts <- Matrix::t(lin_mat)
all_data[["keep"]] <- keep_vec
all_data <- subset(all_data, keep == 1)

# now, assign a cell to only one lineage via tf-idf
lin_mat <- all_data[["lineage"]]@counts
quantile(lin_mat@x)
quantile(colSums(lin_mat))
tfidf_res <- Signac::RunTFIDF(lin_mat)
# count_vec <- sparseMatrixStats::rowSums2(lin_mat)
for(i in 1:ncol(lin_mat)){
  if(i %% floor(ncol(lin_mat)/10) == 0) cat('*')
  
  lin_idx <- .nonzero_col(lin_mat, i)
  if(length(lin_idx) > 1){
    max_idx <- which.max(tfidf_res[lin_idx, i]) # count_vec[lin_idx])
    lin_mat[lin_idx[-max_idx],i] <- 0
  }
}
lin_mat_safe <- lin_mat
quantile(colSums(lin_mat))
quantile(rowSums(lin_mat))
idx <- which(rowSums(lin_mat) != 0)
lin_mat <- lin_mat[idx,]

# now go back and check that each lineage only has one cell-type
lin_mat <- Matrix::t(lin_mat)
idx_list <- lapply(1:ncol(lin_mat), function(j){
  .nonzero_col(lin_mat, j)
})
names(idx_list) <- colnames(lin_mat)
factor_vec <- as.factor(all_data@meta.data$Original_condition)
tabulate_mat <- t(sapply(idx_list, function(idx){
  table(factor_vec[idx])
}))
rownames(tabulate_mat) <- colnames(lin_mat)
naive_idx <- which(colnames(tabulate_mat) == "naive")
tabulate_mat <- tabulate_mat[-which(tabulate_mat[,naive_idx] == 0),]
percentage_mat <- tabulate_mat[,-naive_idx]
percentage_mat <- percentage_mat[rowSums(percentage_mat) != 0,]
percentage_mat <- .mult_vec_mat(1/rowSums(percentage_mat), percentage_mat)

keep_lin_list <- lapply(1:ncol(percentage_mat), function(k){
  rownames(percentage_mat)[which(percentage_mat[,k] > 0.9)]
})
names(keep_lin_list) <- colnames(percentage_mat)
sapply(keep_lin_list, length)

all_lin <- unlist(keep_lin_list)
lin_mat <- lin_mat[,all_lin]
dim(lin_mat)

quantile(colSums(lin_mat))
quantile(rowSums(lin_mat))
keep_vec <- rep(0, nrow(lin_mat))
keep_vec[sparseMatrixStats::rowSums2(lin_mat) > 0] <- 1
table(keep_vec)

all_data[["lineage"]]@data <- Matrix::t(lin_mat)
all_data[["lineage"]]@counts <- Matrix::t(lin_mat)
all_data[["keep"]] <- keep_vec
all_data <- subset(all_data, keep == 1)

#####
lin_mat <- all_data[["lineage"]]@counts
quantile(rowSums(lin_mat))
quantile(colSums(lin_mat))

save(all_data, keep_lin_list,
     file = "../../../../out/kevin/Writeup3d/09302021_sydney_de.RData")
