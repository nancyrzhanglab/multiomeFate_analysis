rm(list=ls())
library(Seurat)
library(multiomeFate)
library(fastTopics)

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup11b/"

load(paste0(out_folder, "Writeup11b_larry-normalized.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- SeuratObject::LayerData(object = seurat_obj, 
                               assay = "RNA", 
                               layer = "counts")
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

save(topic_res, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11b_fasttopics_result.RData"))

#########

topic_mat <- topic_res$L
rownames(topic_mat) <- Seurat::Cells(seurat_obj)
seurat_obj[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                          loadings =  topic_res$F,
                                                          assay = "RNA",
                                                          key =  "fastTopic_")

###########

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(seurat_obj, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11b_larry_fasttopics.RData"))

print("Done! :)")
