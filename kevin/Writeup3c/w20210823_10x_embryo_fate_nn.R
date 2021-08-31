rm(list=ls())
load("../../../../out/kevin/Writeup3c/20210823_10x_embryo_result.RData")

library(Seurat); library(Signac); library(igraph)
source("nn_plot.R")

end_states <- list(16, 6, c(1,2,4), 9)
names(end_states) <- c("Oligo.", "Forebrain", "Cortical1", "Cortical2")
celltype <- mbrain3@meta.data$new_seurat_clusters

#################

tmp <- nn_plot(mbrain3, res$nn_mat, celltype, end_states)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_nndist.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

#############

tmp <- dirnn_plot(mbrain3, res$list_diagnos, celltype, end_states)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_dirnndist.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

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

########################

# try using different distance metric

mat_x <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))
peak_names <- sort(unique(mbrain3[["ATAC"]]@links$peak))
gene_names <- sort(unique(mbrain3[["ATAC"]]@links$gene))
mat_x <- mat_x[,which(colnames(mat_x) %in% peak_names)]
mat_y <- mat_y[,which(colnames(mat_y) %in% gene_names)]

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
                                                                     nn_include_y = F))

tmp <- nn_plot(mbrain3, prep_obj$nn_mat, celltype, end_states)
igraph::diameter(tmp$g, directed = F) #10
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_nndist_euclidean.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")


############################################

# try using both RNA and ATAC for the DE genes

# grab the relevant genes and peaks
mat_x <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))
peak_names <- sort(unique(mbrain3[["ATAC"]]@links$peak))
gene_names <- sort(unique(mbrain3[["ATAC"]]@links$gene))
mat_x <- mat_x[,which(colnames(mat_x) %in% peak_names)]
mat_y <- mat_y[,which(colnames(mat_y) %in% gene_names)]

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
igraph::diameter(tmp$g, directed = F) #10
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_nndist_both.png",
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
igraph::diameter(tmp$g, directed = F) #10
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_de_nndist_both_euclidean.png",
                   tmp$cowplot, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

###########################################
# compute some new dim-reducs

set.seed(10)
tmp <- compute_nndist(mat_x, mat_y, options = list(nn_nn = 5, nn_metric = "euclidean",
                                                   dim_dims_x = 2:rank_x,
                                                   dim_dims_y = 1:rank_y, 
                                                   nn_include_x = T,
                                                   nn_include_y = T))
set.seed(10)
tmp2 <- Seurat::RunUMAP(tmp$all_scores)
tmp2 <- tmp2@cell.embeddings
rownames(tmp2) <- rownames(mbrain3@meta.data)
colnames(tmp2) <- paste0("newwnn_", 1:ncol(tmp2))
mbrain3[["newwnn"]] <- Seurat::CreateDimReducObject(embedding = tmp2, key = "newwnn_")

plot1 <- Seurat::DimPlot(mbrain3, reduction = "newwnn", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nNew Joint (Celltype)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de2_wnn.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")


plot1 <- Seurat::DimPlot(mbrain3, reduction = "newwnn", group.by = "new_seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nNew Joint (Seurat clusters)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de2_wnn_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")

#######

set.seed(10)
tmp <- compute_nndist(mat_x, mat_y, options = list(nn_nn = 5, nn_metric = "euclidean",
                                                   dim_dims_x = 2:rank_x,
                                                   dim_dims_y = 1:rank_y, 
                                                   nn_include_x = T,
                                                   nn_include_y = F))
set.seed(10)
tmp2 <- Seurat::RunUMAP(tmp$all_scores)
tmp2 <- tmp2@cell.embeddings
rownames(tmp2) <- rownames(mbrain3@meta.data)
colnames(tmp2) <- paste0("newwnn_", 1:ncol(tmp2))
mbrain3[["newwnn"]] <- Seurat::CreateDimReducObject(embedding = tmp2, key = "newwnn_")

plot1 <- Seurat::DimPlot(mbrain3, reduction = "newwnn", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nNew ATAC (Celltype)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de2_atac.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")


plot1 <- Seurat::DimPlot(mbrain3, reduction = "newwnn", group.by = "new_seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5, raster = F)
plot1 <- plot1 + ggplot2::ggtitle("Mouse embryo (10x)\nNew ATAC (Seurat clusters)")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_de2_atac_cluster.png", 
                plot1, device = "png", width = 6, height = 5, units = "in")



