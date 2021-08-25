rm(list=ls())
load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed.RData")

# now we want to remove certain cells in "Forebrain GABAergic" or "Cortical or hippocampal glutamatergic"
# NOTE: be sure to check the specific clusters prior to removing these clusters
#   I've noticed that sometimes setting the seed isn't enough to guarantee the same exact cluster-labelings
cell_idx <- which(!mbrain2@meta.data$seurat_clusters %in% c(5, 7, 10, 12, 13, 14))

mat_y <- Matrix::t(mbrain2[["SCT"]]@data[Seurat::VariableFeatures(mbrain2, assay = "SCT"),])[cell_idx,gene_reparam_mat$org_idx]
mat_x <- Matrix::t(mbrain2[["ATAC"]]@data)[cell_idx,peak_reparam_mat$org_idx]

celltype <- mbrain2@meta.data$seurat_clusters[cell_idx]
vec_start <- which(celltype == 3) #glioblast
list_end <- list(which(celltype == 16), #oligodendrocyte
                 which(celltype == 6), #forebrain gabaergic
                 which(celltype %in% c(1,2,4)), #one of cortical
                 which(celltype %in% 9)) #one of cortical
(length(unlist(list_end))+length(vec_start))/length(celltype)

rank_x <- 50
rank_y <- 50
df_x <- data.frame(name = colnames(mat_x))
df_y <- data.frame(name = colnames(mat_y))
p1 <- ncol(mat_x); p2 <- ncol(mat_y)
mat_x <- as.matrix(mat_x)
mat_y <- as.matrix(mat_y)

# form the hash mapping
ht_map <- hash::hash()
for(i in 1:p2){
  ht_map[[as.character(i)]] <- relevant_peaks2[[i]]
}

set.seed(10)
prep_obj <- multiomeFate::chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                                      vec_start = vec_start, list_end = list_end,
                                                      form_method = "average_weighted",
                                                      est_method = "threshold_glmnet",
                                                      cand_method = "nn_any",
                                                      rec_method = "distant_cor",
                                                      ht_map = ht_map,
                                                      options = list(nn_nn = 10, nn_metric = "cosine",
                                                                     dim_dims_x = 2:rank_x,
                                                                     dim_dims_y = 1:rank_y, 
                                                                     est_num_iterations = 4,
                                                                     rec_bool_pred_nn = T,
                                                                     est_cv_choice = "lambda.min",
                                                                     form_bool_include_start = F,
                                                                     form_stepsize = 0.5,
                                                                     form_min_weight = 0,
                                                                     est_verbose = T,
                                                                     rec_verbose = T))

##########

set.seed(10)
res <- multiomeFate::chromatin_potential(prep_obj, verbose = T, bool_oracle = F,
                                         filepath = "../../../../out/kevin/Writeup3c/20210819_10x_embryo_result_tmp.RData")


save.image("../../../../out/kevin/Writeup3c/20210819_10x_embryo_result.RData")

