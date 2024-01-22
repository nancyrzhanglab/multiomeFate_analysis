rm(list=ls())
library(Seurat)
library(SAVER)

load("../../../../out/kevin/analysis_pipeline/step3_basic-eda.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
saver_res <- SAVER::saver(x = mat, ncores = 4)

#######

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

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/analysis_pipeline/step4_saver.RData")

