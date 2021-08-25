rm(list=ls())
load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed_de.RData")

library(Seurat); library(Signac)

# grab the relevant genes and peaks
mat_x <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))
peak_names <- sort(unique(mbrain3[["ATAC"]]@links$peak))
gene_names <- sort(unique(mbrain3[["ATAC"]]@links$gene))
mat_x <- mat_x[,which(colnames(mat_x) %in% peak_names)]
mat_y <- mat_y[,which(colnames(mat_y) %in% gene_names)]

# [[note to self 1: can a peak be associated with multiple genes?]]
# [[note to self 2: check to make sure the genes are indeed DE]]

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

############

celltype <- mbrain3@meta.data$new_seurat_clusters
vec_start <- which(celltype == 3) #glioblast
list_end <- list(which(celltype == 16), #oligodendrocyte
                 which(celltype == 6), #forebrain gabaergic
                 which(celltype %in% c(1,2,4)), #one of cortical
                 which(celltype %in% 9)) #one of cortical

############

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
                                                                     nn_include_y = F,
                                                                     est_num_iterations = 4,
                                                                     rec_bool_pred_nn = T,
                                                                     est_cv_choice = "lambda.min",
                                                                     form_bool_include_start = F,
                                                                     form_stepsize = 1,
                                                                     est_verbose = T,
                                                                     rec_verbose = T))

##########

time_start <- Sys.time()
set.seed(10)
res <- multiomeFate::chromatin_potential(prep_obj, verbose = T, bool_oracle = F,
                                         filepath = "../../../../out/kevin/Writeup3c/20210823_10x_embryo_result_tmp.RData")
time_end <- Sys.time()

save.image("../../../../out/kevin/Writeup3c/20210823_10x_embryo_result.RData")

