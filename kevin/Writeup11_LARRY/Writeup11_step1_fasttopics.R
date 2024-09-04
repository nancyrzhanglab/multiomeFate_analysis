rm(list=ls())
library(Seurat)
library(fastTopics)
load("~/nzhanglab/data/LARRY_pyro-velocity/Writeup11_larry_multilineage.RData")

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup11/"

mat <- SeuratObject::LayerData(object = seurat_obj, 
                               assay = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(seurat_obj, assay = "RNA"))
mat <- Matrix::t(mat)

K <- 15
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

save(topic_res, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11_larry_fasttopics.RData"))
