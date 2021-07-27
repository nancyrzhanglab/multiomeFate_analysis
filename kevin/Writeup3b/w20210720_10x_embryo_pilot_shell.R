
rm(list=ls()); gc(T)
load("../../../../out/kevin/Writeup3b/10x_embryo.RData")

vec_start <- which(celltype == "Radial glia")
list_end <- list(which(celltype == "Oligodendrocyte"),
                 which(celltype == "Forebrain GABAergic"),
                 which(celltype == "Cortical or hippocampal glutamatergic"))
(length(vec_start) + length(unlist(list_end)))/nrow(mat_x)

rank_x <- 50
rank_y <- 30
df_x <- data.frame(name = colnames(mat_x))
df_y <- data.frame(name = colnames(mat_y))
mat_x <- as.matrix(mat_x)
mat_y <- as.matrix(mat_y)

set.seed(10)
prep_obj <- multiomeFate::chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                                      vec_start, list_end,
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
                                                                     form_stepsize = 0.5,
                                                                     form_min_weight = 0,
                                                                     est_verbose = T,
                                                                     rec_verbose = T))

set.seed(10)
res <- multiomeFate::chromatin_potential(prep_obj, verbose = T, bool_oracle = F,
                                         filepath = "../../../../out/kevin/Writeup3b/10x_embryo_result_tmp.RData")

save.image("../../../../out/kevin/Writeup3b/10x_embryo_result2.RData")