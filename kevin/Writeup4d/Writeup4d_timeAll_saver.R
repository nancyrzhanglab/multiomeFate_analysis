rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_peakmerging.RData")
library(Seurat)
library(SAVER)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
print(dim(mat))
saver_res <- SAVER::saver(x = mat, ncores = 4)
save(saver_res, all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver.RData")

print("Computing UMAP")
tmp <- saver_res$estimate
tmp <- pmin(tmp, 10)
all_data[["Saver"]] <- Seurat::CreateAssayObject(counts = tmp)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data <- Seurat::ScaleData(all_data)
all_data[["Saver"]]@var.features <- rownames(tmp)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE,
                           reduction.name = "saverpca") 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 1:50,
                            reduction = "saverpca",
                            reduction.name = "saverumap")
save(saver_res, all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver.RData")


