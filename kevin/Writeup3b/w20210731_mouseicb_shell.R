rm(list=ls())
load("../../../../out/kevin/Writeup3b/mouseicb_fate_prep.RData")
load("../../../../out/kevin/Writeup3b/20210731_mouseicb_preprocess2.RData")

Seurat::DefaultAssay(myeloid2) <- "RNA"
myeloid2 <- Seurat::NormalizeData(myeloid2)
mat_y <- as.matrix(Matrix::t(myeloid2[["RNA"]]@data))
mat_x <- Matrix::t(myeloid2[["ATAC"]]@data)

celltype <- myeloid2@meta.data$celltype2
list_end <- list(which(celltype == "Cluster2"),
                 which(celltype == "Cluster3"),
                 which(celltype == "Cluster4"),
                 which(celltype == "Cluster5"),
                 which(celltype == "Cluster6"))

rank_x <- 50
rank_y <- 30
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
                                                      vec_start = NA, list_end,
                                                      form_method = "average_weighted",
                                                      est_method = "threshold_glmnet",
                                                      cand_method = "nn_freq",
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

set.seed(10)
res <- multiomeFate::chromatin_potential(prep_obj, verbose = T, bool_oracle = F,
                                         filepath = "../../../../out/kevin/Writeup3b/20210731_mouseicb_result_tmp.RData")


save.image("../../../../out/kevin/Writeup3b/20210731_mouseicb_result.RData")