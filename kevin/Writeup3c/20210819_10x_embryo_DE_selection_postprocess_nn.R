rm(list=ls())
load("../../../../out/kevin/Writeup3c/20210823_10x_embryo_result.RData")

library(Seurat); library(Signac)

n <- nrow(mat_x)
mat <- matrix(0, n, n)
for(i in 1:n){
  mat[i,res$nn_mat[i,]] <- 1
}

library(igraph)
g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected")
g <- igraph::simplify(g)

igraph::radius(g) # radius of 7...

celltype <- mbrain3@meta.data$new_seurat_clusters
oligo_idx <- which(celltype == 16)
dist_oligo <- igraph::distances(g, to = oligo_idx)
tmp1 <- apply(dist_oligo, 1, mean)
tmp2 <- apply(dist_oligo, 1, min)
tmp3 <- apply(dist_oligo, 1, max)
dist_oligo_summary <- cbind(tmp1, tmp2, tmp3)
colnames(dist_oligo_summary) <- c("oligo_mean", "oligo_min", "oligo_max")
rownames(dist_oligo_summary) <- rownames(mbrain3@meta.data)

########
apply(dist_oligo_summary[oligo_idx,], 2, mean)
glio_idx <- which(celltype == 3) 
apply(dist_oligo_summary[glio_idx,], 2, mean)
forebrain_idx <- which(celltype == 6) 
apply(dist_oligo_summary[forebrain_idx,], 2, mean)
cortical1_idx <- which(celltype %in% c(1,2,4)) 
apply(dist_oligo_summary[cortical1_idx,], 2, mean)


