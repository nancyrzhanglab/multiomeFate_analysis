rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)
library(fastTopics)

load("../../../../out/kevin/Writeup5a/Writeup5a_tcca_RNA-geneActivity.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)
note <- "Processed on Nov-5-2022. The fasttopics need to be replaced, as they were computed when we still used the Test2 and Test3 datasets (see Writeup4d_timeAll_fasttopics_CIS.R). Also, the tilted-CCA outputs here were previously computed when gene activities were only +2000bp up, 0bp"

print("Seurat preprocess of gene activity")
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["geneActivity"]] <- NULL

gene_activities <- Signac::GeneActivity(all_data)
all_data[["geneActivity"]] <- CreateAssayObject(counts = gene_activities,
                                                extend.upstream = 1000,
                                                extend.downstream = 1000)
Seurat::DefaultAssay(all_data) <- "geneActivity"
all_data[["geneActivity"]]@var.features <- intersect(rownames(gene_activities), all_data[["RNA"]]@var.features)
all_data <- Seurat::NormalizeData(
  object = all_data,
  assay = "geneActivity",
  normalization.method = "LogNormalize"
)
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE,
                           reduction.name = "activityPCA") 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 2:50,
                            reduction = "activityPCA", 
                            reduction.name = "activity.umap")


print("Finished preprocessing data")
save(all_data, date_of_run, session_info, note,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

