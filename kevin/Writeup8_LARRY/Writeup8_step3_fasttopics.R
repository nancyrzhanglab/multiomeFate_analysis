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
sum_vec <- Matrix::rowSums(mat)
if(any(sum_vec <= 0.3)){
  mat <- mat[sum_vec > 0.3, ]
}
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

#########

topic_mat <- matrix(NA, nrow = ncol(seurat_object), ncol = ncol(topic_res$L))
rownames(topic_mat) <- colnames(seurat_object)
topic_res$L <- topic_res$L[rownames(topic_res$L) %in% colnames(seurat_object),]

for(i in 1:nrow(topic_res$L)){
  topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
}
colnames(topic_mat) <- paste0("fastTopic_", 1:ncol(topic_mat))
colnames(topic_res$F) <- paste0("fastTopic_", 1:ncol(topic_mat))

seurat_object[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                             loadings =  topic_res$F,
                                                             assay = "RNA",
                                                             key =  "fastTopic_")

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(seurat_object, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step3_fasttopics.RData")

print("Done! :)")
