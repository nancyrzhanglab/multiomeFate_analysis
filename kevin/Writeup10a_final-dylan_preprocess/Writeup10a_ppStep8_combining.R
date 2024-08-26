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
  full_umap_res <- matrix(NA, nrow = length(cell_names), ncol = 2)
  rownames(topic_mat) <- cell_names
  rownames(full_umap_res) <- cell_names
  
  set.seed(10)
  umap_res <- Seurat::RunUMAP(topic_res$L)@cell.embeddings
  
  for(i in 1:nrow(topic_res$L)){
    topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
    full_umap_res[rownames(umap_res)[i],] <- umap_res[i,]
  }
  colnames(topic_mat) <- paste0("fastTopic", treatment, "_", 1:K)
  colnames(topic_res$F) <- paste0("fastTopic", treatment, "_", 1:K)
  colnames(full_umap_res) <- paste0("ft", treatment, "umap_", 1:2)
  
  all_data[[paste0("fasttopic.", treatment)]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                                              loadings =  topic_res$F,
                                                                              assay = "RNA")
  all_data[[paste0("ft.", treatment, ".umap")]] <- Seurat::CreateDimReducObject(embeddings = full_umap_res,
                                                                                assay = "RNA")
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
                            reduction = "saver.pca",
                            reduction.name = "saver.umap")

#################

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

print("Saving")
save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep8_combining.RData"))

print("Done! :)")

