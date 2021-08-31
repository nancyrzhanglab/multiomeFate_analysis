rm(list=ls())
load("../../../../out/kevin/Writeup3c/20210823_10x_embryo_result.RData")

library(Seurat); library(Signac); library(igraph)
source("nn_plot.R")

end_states <- list(16, 6, c(1,2,4), 9)
names(end_states) <- c("Oligo.", "Forebrain", "Cortical1", "Cortical2")
celltype <- mbrain3@meta.data$new_seurat_clusters

############################################
# all the things we need to run the multiome_prep

# create the hash map
p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
ht_map <- hash::hash()
for(i in 1:p2){
  gene_name <- colnames(mat_y)[i]
  idx <- which(mbrain3[["ATAC"]]@links$gene == gene_name)
  peak_names <- mbrain3[["ATAC"]]@links$peak[idx]
  ht_map[[as.character(i)]] <- which(colnames(mat_x) %in% peak_names)
}

rank_x <- 50
rank_y <- 30
df_x <- data.frame(name = colnames(mat_x))
df_y <- data.frame(name = colnames(mat_y))

celltype <- mbrain3@meta.data$new_seurat_clusters
vec_start <- which(celltype == 15) #radial
list_end <- list(which(celltype == 16), #oligodendrocyte
                 which(celltype == 6), #forebrain gabaergic
                 which(celltype %in% c(1,2,4)), #one of cortical
                 which(celltype %in% 9)) #one of cortical

#############################

mat_x <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))
peak_names <- sort(unique(mbrain3[["ATAC"]]@links$peak))
gene_names <- sort(unique(mbrain3[["ATAC"]]@links$gene))
mat_x <- mat_x[,which(colnames(mat_x) %in% peak_names)]
mat_y <- mat_y[,which(colnames(mat_y) %in% gene_names)]

set.seed(10)

# determine the densities
nn_vec <- c(5, 10, 20, 50)
tmp <- density_plot(mbrain3, mat_x, mat_y,
                    nn_vec, include_x = T, include_y = F)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_density.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")


##############

nn_vec <- c(5, 10, 20, 50)
tmp <- density_plot(mbrain3, mat_x, mat_y,
                    nn_vec, include_x = T, include_y = T, nn_metric = "euclidean")
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_density_both_euclidean.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

#####################################

set.seed(10)
tmp <- compute_nndist(mat_x, mat_y, options = list(nn_nn = 5, nn_metric = "cosine",
                                                   dim_dims_x = 2:rank_x,
                                                   dim_dims_y = 1:rank_y, 
                                                   nn_include_x = T,
                                                   nn_include_y = F))
dim(tmp$all_scores); range(tmp$all_scores)
idx <- which(celltype == "16"); quantile(tmp$vec[idx])

set.seed(10)
tmp <- compute_nndist(mat_x, mat_y, options = list(nn_nn = 5, nn_metric = "euclidean",
                                                   dim_dims_x = 2:rank_x,
                                                   dim_dims_y = 1:rank_y, 
                                                   nn_include_x = T,
                                                   nn_include_y = T))
dim(tmp$all_scores); range(tmp$all_scores)
idx <- which(celltype == "16"); quantile(tmp$vec[idx])
idx <- which(celltype == "9"); quantile(tmp$vec[idx])




