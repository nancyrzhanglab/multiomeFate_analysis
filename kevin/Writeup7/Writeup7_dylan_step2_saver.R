rm(list=ls())
library(Seurat)
library(SAVER)

load("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step1_basic-eda.RData")
set.seed(10)

Seurat::DefaultAssay(all_data) <- "RNA"
mat <- Seurat::GetAssayData(object = all_data, assay = "RNA", slot = "counts")
mat <- mat[Seurat::VariableFeatures(all_data),]
print(dim(mat))
print(mat[1:5,1:5])

print("Running SAVER")
saver_res <- SAVER::saver(x = mat, ncores = 4)

save(saver_res,
     file = "~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step2_saver.RData")

#######

print("Postprocessing SAVER results")
tmp <- saver_res$estimate
tmp <- tmp[,colnames(all_data)]
tmp <- pmin(tmp, 10)
all_data[["Saver"]] <- Seurat::CreateAssayObject(data = tmp)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data <- Seurat::ScaleData(all_data)
Seurat::VariableFeatures(all_data) <- rownames(tmp)
all_data <- Seurat::RunPCA(all_data, 
                           assay = "Saver",
                           verbose = FALSE,
                           reduction.name = "saverpca") 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            assay = "Saver",
                            dims = 1:50,
                            reduction = "saverpca",
                            reduction.name = "saverumap")

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(all_data, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step2_saver.RData")

print("Done! :)")
