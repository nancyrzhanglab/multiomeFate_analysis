rm(list=ls())
load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed_de.RData")

library(Seurat); library(Signac)

Seurat::DefaultAssay(mbrain3) <- "ATAC"
mbrain3 <- Seurat::FindNeighbors(
  object = mbrain3,
  reduction = 'lsi',
  dims = 2:50
)

set.seed(10)
mbrain3 <- Seurat::FindClusters(
  object = mbrain3,
  algorithm = 3,
  resolution = 10,
  verbose = T
)

############################

reduc_vec <- c("umap.rna", "umap.atac", "wnn.umap")
main_vec <- c("RNA", "ATAC", "WNN")
file_vec <- c("rna", "atac", "wnn")

uniq_clust <- sort(unique(as.character(mbrain3@meta.data$ATAC_snn_res.10)))
col_vec <- scales::hue_pal()(length(uniq_clust))

for(i in 1:3){
  plot1 <- Seurat::DimPlot(mbrain3, reduction = reduc_vec[i], group.by = "ATAC_snn_res.10", label = TRUE,
                           repel = TRUE, label.size = 2.5, raster = F)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x)\n", main_vec[i], " (Meta-cells clustering by Seurat)"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_", file_vec[i], "_metacell_clustering.png"), 
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  png(paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_", file_vec[i], "_metacell.png"),
      height = 1500, width = 1500, res = 300, units = "px")
  graphics::plot(NA, 
                 xlim = range(mbrain3[[reduc_vec[i]]]@cell.embeddings[,1]),
                 ylim = range(mbrain3[[reduc_vec[i]]]@cell.embeddings[,2]),
                 xlab = paste0(file_vec[i], "UMAP_1"),
                 ylab = paste0(file_vec[i], "UMAP_2"),
                 main = paste0("Mouse embryo (10x)\n", main_vec[i], " (Meta-cells by Seurat)"))
  for(k in 1:length(uniq_clust)){
    clust <- uniq_clust[k]
    idx <- which(as.character(mbrain3@meta.data$ATAC_snn_res.10) == clust)
    coord <- apply(mbrain3[[reduc_vec[i]]]@cell.embeddings[idx,], 2, median)
    points(coord[1], coord[2], 
           pch = 16, 
           cex = 2, 
           col = col_vec[k])
  }
  graphics.off()
}

#####################################

uniq_clust <- sort(unique(as.character(mbrain3@meta.data$ATAC_snn_res.10)))
lsi_mat <- mbrain3[["lsi"]]@cell.embeddings
clust_mat <- t(sapply(uniq_clust, function(clust){
  idx <- which(as.character(mbrain3@meta.data$ATAC_snn_res.10) == clust)
  apply(lsi_mat[idx,,drop = F], 2, median)
}))
rownames(clust_mat) <- uniq_clust

set.seed(10)
nn_res <- knn.covertree::find_knn(clust_mat, 
                                  k = 10,
                                  distance = "cosine")

n <- nrow(clust_mat)
adj_mat <- matrix(0, n, n)
for(i in 1:n){
  adj_mat[i, nn_res$index[i,]] <- 1
}
diag(adj_mat) <- 0
adj_mat <- adj_mat*t(adj_mat)
idx <- which(rowSums(adj_mat) == 0)
for(i in idx){
  adj_mat[i,nn_res$nn.idx[i,2]] <- 1
}
adj_mat <- adj_mat + t(adj_mat)
adj_mat[adj_mat > 0] <- 1
diag(adj_mat) <- 0

for(i in 1:3){
  median_coord <- t(sapply(uniq_clust, function(clust){
    idx <- which(as.character(mbrain3@meta.data$ATAC_snn_res.10) == clust)
    apply(mbrain3[[reduc_vec[i]]]@cell.embeddings[idx,,drop = F], 2, median)
  }))
  
  png(paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_", file_vec[i], "_metacell_knn.png"),
      height = 1500, width = 1500, res = 300, units = "px")
  graphics::plot(NA, 
                 xlim = range(mbrain3[[reduc_vec[i]]]@cell.embeddings[,1]),
                 ylim = range(mbrain3[[reduc_vec[i]]]@cell.embeddings[,2]),
                 xlab = paste0(file_vec[i], "UMAP_1"),
                 ylab = paste0(file_vec[i], "UMAP_2"),
                 main = paste0("Mouse embryo (10x)\n", main_vec[i], " (Meta-cells by Seurat)"))
  for(j in 1:nrow(clust_mat)){
    idx <- which(adj_mat[j,] != 0)
    for(j2 in idx){
      lines(median_coord[c(j,j2),])
    }
  }
  
  for(k in 1:length(uniq_clust)){
    clust <- uniq_clust[k]
    idx <- which(as.character(mbrain3@meta.data$ATAC_snn_res.10) == clust)
    coord <- apply(mbrain3[[reduc_vec[i]]]@cell.embeddings[idx,], 2, median)
    points(coord[1], coord[2], 
           pch = 16, 
           cex = 2, 
           col = col_vec[k])
  }
  graphics.off()
}

################

# try the resistance distance
# see https://www.stat.berkeley.edu/~mmahoney/s15-stat260-cs294/Lectures/lecture16-17mar15.pdf
laplacian <- diag(colSums(adj_mat)) - adj_mat
lap_inv <- MASS::ginv(laplacian)
start <- "4"
end_vec <- c("46", "11", "1") # oligo, forebrain, cortical
start_idx <- which(uniq_clust == start)
n <- length(uniq_clust)
sapply(end_vec, function(end){
  end_idx <- which(uniq_clust == end)
  vec <- rep(0, n)
  vec[start_idx] <- 1
  vec[end_idx] <- -1
  t(vec) %*% lap_inv %*% vec
})

# try the diffusion distance, see 
# https://www.stat.berkeley.edu/~mmahoney/s15-stat260-cs294/Lectures/lecture15-12mar15.pdf
# https://academic.oup.com/bioinformatics/article/31/18/2989/241305
# https://arxiv.org/pdf/1406.0013.pdf
# https://mathworld.wolfram.com/LeftEigenvector.html
# left and right eigenvalues
n <- nrow(P)
P <- adj_mat
tmp <- matrixStats::rowSums2(P)
P <- diag(1/tmp) %*% P %*% diag(1/tmp)
P <- diag(1/matrixStats::rowSums2(P)) %*% P
right_eigen <- eigen(P, symmetric = F)
left_eigen <- eigen(t(P), symmetric = F)
if(any(left_eigen$vectors[,1] < 0)) left_eigen$vectors[,1] <- -left_eigen$vectors[,1]
eigenvalues <- left_eigen$values 
# check
sum(abs(P%*%right_eigen$vectors[,2] - right_eigen$values[2]*right_eigen$vectors[,2]))
sum(abs(left_eigen$vectors[,2]%*%P - left_eigen$values[2]*left_eigen$vectors[,2]))
sum(abs(right_eigen$vectors[,1]+rep(1/sqrt(n))))
sum(abs(right_eigen$values - left_eigen$values))
# normalize
sum(left_eigen$vectors[,2]^2)
left_vector <- left_eigen$vectors
for(i in 2:ncol(left_vector)){
  left_vector[,i] <- left_vector[,i]*sqrt(left_vector[,1])
}
right_vector <- right_eigen$vectors
for(i in 1:ncol(right_vector)){
  right_vector[,i] <- right_vector[,i]/sqrt(left_vector[,1])
}
sum(left_vector[,2]^2/left_vector[,1])
sum(right_vector[,2]^2*left_vector[,1])

diffusion_distance <- function(eigenvalues, right_vector, 
                               idx1, idx2, 
                               time_vec = 1:40){
  sqrt(sum(sapply(time_vec, function(x){
    eigenvalues^(2*x)*sum((right_vector[idx1,-1] - right_vector[idx2,-1])^2)
  })))
}


cluster_vec <- c("4", # radial
                 "46", # oligo
                 "11", # forebrain
                 "1", "8", #cortical1, cortical2
                 "31", "5", # neuroblast1, neuroblast2
                 "27") #glio
n <- length(uniq_clust)
dist_mat <- matrix(0, length(cluster_vec), length(cluster_vec))
colnames(dist_mat) <- rownames(dist_mat) <- cluster_vec
for(i in 1:(length(cluster_vec)-1)){
  idx1 <- which(uniq_clust == cluster_vec[i])
  for(j in (i+1):length(cluster_vec)){
    idx2 <-  which(uniq_clust == cluster_vec[j])
    dist_mat[i,j] <- diffusion_distance(eigenvalues, right_vector,
                       idx1, idx2)
  }
}
round(dist_mat,2)
