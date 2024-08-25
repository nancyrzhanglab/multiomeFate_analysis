rm(list=ls())
library(Seurat)
library(Signac)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc.RData"))
treatment_vec <- c("CIS", "COCL2", "DABTRAM")

# process  fasttopics

print("Working on fasttopics")
for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  load(paste0(out_folder, "Writeup10a_ppStep7a_fasttopic_", treatment, ".RData"))
  
  cell_names <- Seurat::Cells(all_data)
  K <- ncol(topic_res$L)
  topic_mat <- matrix(NA, nrow = length(cell_names), ncol = K)
  rownames(topic_mat) <- cell_names
  topic_res$L <- topic_res$L[rownames(topic_res$L) %in% cell_names,]
  
  for(i in 1:nrow(topic_res$L)){
    topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
  }
  colnames(topic_mat) <- paste0("fastTopic", treatment, "_", 1:K)
  colnames(topic_res$F) <- paste0("fastTopic", treatment, "_", 1:K)
  
  all_data[[paste0("fasttopic_", treatment)]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                                              loadings =  topic_res$F,
                                                                              assay = "RNA",
                                                                              key =  paste0("fastTopic", treatment, "_"))
  
  print("Computing UMAP")
  set.seed(10)
  all_data <- Seurat::RunUMAP(all_data, 
                              dims = 1:K,
                              reduction = paste0("fasttopic_", treatment),
                              reduction.name = paste0("ft.", treatment, ".umap"))
                                                                              
}

#############################

# process SAVER
print("Working on SAVER")
load(paste0(out_folder, "Writeup10a_ppStep7b_saver.RData"))

tmp <- saver_res$estimate
tmp <- tmp[,Seurat::Cells(all_data)]
tmp <- pmin(tmp, 10)
all_data[["Saver"]] <- Seurat::CreateAssayObject(counts = tmp)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data <- Seurat::ScaleData(all_data)
Seurat::VariableFeatures(all_data[["Saver"]]) <- rownames(tmp)
all_data <- Seurat::RunPCA(all_data, 
                           verbose = FALSE,
                           reduction.name = "saver.pca")
print("Computing UMAP")
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            dims = 1:50,
                            reduction = "saverpca",
                            reduction.name = "saver.umap")

#################

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

print("Saving")
save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep8_combining.RData"))

print("Done! :)")

