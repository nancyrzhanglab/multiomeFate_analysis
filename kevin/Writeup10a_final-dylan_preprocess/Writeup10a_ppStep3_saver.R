rm(list=ls())
library(Seurat)
library(Signac)
library(SAVER)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
print(dim(mat))
saver_res <- SAVER::saver(x = mat, ncores = 6)
save(saver_res, all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_saver.RData")

print("Computing UMAP")
keep_vec <- rep(0, ncol(all_data))
keep_vec[which(colnames(all_data) %in% colnames(saver_res$estimate))] <- 1
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == 1)

tmp <- saver_res$estimate
tmp <- tmp[,colnames(all_data)]
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
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_saver.RData")
