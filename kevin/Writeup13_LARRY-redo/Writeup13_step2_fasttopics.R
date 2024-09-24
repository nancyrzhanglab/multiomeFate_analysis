rm(list=ls())
library(Seurat)
library(multiomeFate)
load("~/project/Multiome_fate/out/kevin/Writeup13/Writeup13_larry-dataset.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(seurat_object) <- "RNA"
mat <- SeuratObject::LayerData(object = seurat_object, 
                               assay = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(seurat_object))
sum_vec <- Matrix::colSums(mat)
mat <- multiomeFate:::.mult_mat_vec(mat,1/sum_vec)
mat@x <- round(1e4*mat@x)
sum_vec <- Matrix::colSums(mat)

sum_vec <- Matrix::rowSums(mat)
if(any(sum_vec <= 1)){
  mat <- mat[sum_vec > 1, ]
}
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

topic_mat <- topic_res$L
rownames(topic_mat) <- Seurat::Cells(seurat_object)
seurat_object[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                             loadings =  topic_res$F,
                                                             assay = "RNA",
                                                             key =  "fastTopic_")

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(seurat_object, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup13/Writeup13_larry-dataset_step2_fasttopics.RData")

print("Done! :)")
