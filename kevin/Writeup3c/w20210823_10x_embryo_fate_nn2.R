rm(list=ls())
load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed.RData")
load("../../../../out/kevin/Writeup3c/20210823_10x_embryo_result.RData")

library(Seurat); library(Signac); library(igraph)
source("nn_plot.R")

end_states <- list(16, 6, c(1,2,4), 9)
names(end_states) <- c("Oligo.", "Forebrain", "Cortical1", "Cortical2")
celltype <- mbrain3@meta.data$new_seurat_clusters

mat_y <- Matrix::t(mbrain3[["SCT"]]@data[Seurat::VariableFeatures(mbrain3, assay = "SCT"),])[,gene_reparam_mat$org_idx]
mat_x <- Matrix::t(mbrain3[["ATAC"]]@data)[,peak_reparam_mat$org_idx]
mat_x <- as.matrix(mat_x); mat_y <- as.matrix(mat_y)

##########################

p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
ht_map <- hash::hash()
for(i in 1:p2){
  ht_map[[as.character(i)]] <- relevant_peaks2[[i]]
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

##########################

set.seed(10)
prep_obj <- multiomeFate::chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                                      vec_start = vec_start, list_end = list_end,
                                                      form_method = "average_weighted",
                                                      est_method = "threshold_glmnet",
                                                      cand_method = "nn_any",
                                                      rec_method = "distant_cor",
                                                      ht_map = ht_map,
                                                      options = list(nn_nn = 5, nn_metric = "cosine",
                                                                     dim_dims_x = 2:rank_x,
                                                                     dim_dims_y = 1:rank_y, 
                                                                     nn_include_x = T,
                                                                     nn_include_y = F))
tmp <- nn_plot(mbrain3, prep_obj$nn_mat, celltype, end_states)
igraph::diameter(tmp$g, directed = F) #10
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_nndist.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")


set.seed(10)
prep_obj <- multiomeFate::chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                                      vec_start = vec_start, list_end = list_end,
                                                      form_method = "average_weighted",
                                                      est_method = "threshold_glmnet",
                                                      cand_method = "nn_any",
                                                      rec_method = "distant_cor",
                                                      ht_map = ht_map,
                                                      options = list(nn_nn = 5, nn_metric = "euclidean",
                                                                     dim_dims_x = 2:rank_x,
                                                                     dim_dims_y = 1:rank_y, 
                                                                     nn_include_x = T,
                                                                     nn_include_y = T))
tmp <- nn_plot(mbrain3, prep_obj$nn_mat, celltype, end_states)
igraph::diameter(tmp$g, directed = F) #13
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_nndist_both_euclidean.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")


set.seed(10)
prep_obj <- multiomeFate::chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                                      vec_start = vec_start, list_end = list_end,
                                                      form_method = "average_weighted",
                                                      est_method = "threshold_glmnet",
                                                      cand_method = "nn_any",
                                                      rec_method = "distant_cor",
                                                      ht_map = ht_map,
                                                      options = list(nn_nn = 5, nn_metric = "cosine",
                                                                     dim_dims_x = 2:rank_x,
                                                                     dim_dims_y = 1:rank_y, 
                                                                     nn_include_x = T,
                                                                     nn_include_y = T))
tmp <- nn_plot(mbrain3, prep_obj$nn_mat, celltype, end_states)
igraph::diameter(tmp$g, directed = F) #17
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_nndist_both.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")




