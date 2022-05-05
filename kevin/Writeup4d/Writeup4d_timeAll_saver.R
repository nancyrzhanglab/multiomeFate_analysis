rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_peakmerging.RData")
library(Seurat)
library(Signac)
library(SAVER)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
print(dim(mat))
saver_res <- SAVER::saver(x = mat, ncores = 6)
save(saver_res, all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver.RData")

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
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver.RData")

###############

# visualizations of one treatment are at a time
treatment_vec <- c("CIS", "COCL2", "DABTRAM")

for(treatment in treatment_vec){
  all_data_subset <- all_data
  Seurat::DefaultAssay(all_data_subset) <- "Saver"
  all_data_subset[["ATAC"]] <- NULL
  
  keep_vec <- rep(0, ncol(all_data_subset))
  keep_vec[which(all_data_subset$dataset %in% c("day0", "test2", "test3", 
                                                paste0("day10_", treatment),
                                                paste0("week5_", treatment)))] <- 1
  all_data_subset$keep <- keep_vec
  all_data_subset <- subset(all_data_subset, keep == 1)
  
  all_data_subset[["saverpca"]] <- NULL
  all_data_subset[["saverumap"]] <- NULL
  
  all_data_subset <- Seurat::RunPCA(all_data_subset, verbose = FALSE,
                                    reduction.name = "saverpca") 
  set.seed(10)
  all_data_subset <- Seurat::RunUMAP(all_data_subset, dims = 1:50,
                                     reduction = "saverpca",
                                     reduction.name = "saverumap")
  
  save(all_data_subset, date_of_run, session_info,
       file = paste0("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver_onlyRNA_", treatment, "-subset.RData"))
}



