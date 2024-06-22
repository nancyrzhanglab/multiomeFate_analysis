rm(list=ls())
library(Seurat)
library(multiomeFate)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step2_lineage-plotting.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- Seurat::GetAssayData(object = seurat_object, 
                            assay = "RNA", 
                            layer = "counts")
mat <- mat[Seurat::VariableFeatures(seurat_object),]
sum_vec <- Matrix::colSums(mat)
mat <- multiomeFate:::.mult_mat_vec(mat,1/sum_vec)
mat@x <- round(1e4*mat@x)
sum_vec <- Matrix::colSums(mat)
# quantile(sum_vec)
# quantile(mat@x)

sum_vec <- Matrix::rowSums(mat)
if(any(sum_vec <= 1)){
  mat <- mat[sum_vec > 1, ]
}
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

save(topic_res, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup9c/Writeup9c_larry-dataset_step3_fasttopics_tmp.RData")

#########

topic_mat <- topic_res$L
seurat_object[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                             loadings =  topic_res$F,
                                                             assay = "RNA",
                                                             key =  "fastTopic_")

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(seurat_object, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup9c/Writeup9c_larry-dataset_step3_fasttopics.RData")

print("Done! :)")
