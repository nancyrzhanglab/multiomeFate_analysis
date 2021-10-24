rm(list=ls())
library(Seurat)
load("../../../../out/kevin/Writeup3e/10192021_sydney_preprocess.RData")
source("../Writeup3d/funcs.R")
source("select_cells.R")

###############

Seurat::Idents(all_data) <- "Original_condition"
table(Seurat::Idents(all_data))

tabulate_mat <- .tabulate_lineages(all_data)
tabulate_mat[1:15,]
quantile(tabulate_mat[,"naive"], probs = seq(0,1,length.out=11))
max_val <- sapply(1:nrow(tabulate_mat), function(i){max(tabulate_mat[i,-5])})
tabulate_mat <- tabulate_mat[order(max_val, decreasing = T),]
tabulate_mat[1:15,]
tabulate_mat[nrow(tabulate_mat):1,][1:15,]
tabulate_mat[order(tabulate_mat[,"tram"], decreasing = T),][1:15,]

###################

percentage_mat <- tabulate_mat
percentage_mat <- percentage_mat[rowSums(percentage_mat) != 0,]
percentage_mat <- .mult_vec_mat(1/rowSums(percentage_mat), percentage_mat)

# cluster the percentages
set.seed(10)
eps <- .l2norm(rep(1/6,6) - c(rep(1/5,5),0))/10
dbscan_clust <- dbscan::dbscan(percentage_mat, 
                               eps = 0.1,
                               minPts = 4)
table(dbscan_clust$cluster)

cell_per_cluster <- sapply(sort(unique(dbscan_clust$cluster)), function(clust){
  lineage_idx <- which(dbscan_clust$cluster == clust)
  tmp_mat <- all_data[["lineage"]]@data[rownames(all_data[["lineage"]]@data) %in% rownames(percentage_mat)[lineage_idx],]
  length(which(sparseMatrixStats::colSums2(tmp_mat) > 0))
})
naive_per_cluster <- sapply(sort(unique(dbscan_clust$cluster)), function(clust){
  lineage_idx <- which(dbscan_clust$cluster == clust)
  tmp_mat <- all_data[["lineage"]]@data[rownames(all_data[["lineage"]]@data) %in% rownames(percentage_mat)[lineage_idx],]
  tmp_mat <- tmp_mat[,which(all_data$Original_condition == "naive")]
  length(which(sparseMatrixStats::colSums2(tmp_mat) > 0))
})

zz <- rbind(table(dbscan_clust$cluster), cell_per_cluster, naive_per_cluster)
rownames(zz) <- c("Number-of-lineages", "Number-of-cells", "Number-of-naive-cells")
zz

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

########################


tabulate_mat2 <- tabulate_mat
tabulate_mat2 <- log1p(tabulate_mat2)
tabulate_mat2 <- as.data.frame(tabulate_mat2)
plot1 <- GGally::ggpairs(tabulate_mat2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10192021_sydney_pairs_log.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")


percentage_mat <- tabulate_mat
percentage_mat <- percentage_mat[rowSums(percentage_mat) != 0,]
percentage_mat <- .mult_vec_mat(1/rowSums(percentage_mat), percentage_mat)
percentage_mat <- as.data.frame(percentage_mat)
plot1 <- GGally::ggpairs(percentage_mat)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/10192021_sydney_pairs_percentage.png"),
                plot1, device = "png", width = 7, height = 7, units = "in")

