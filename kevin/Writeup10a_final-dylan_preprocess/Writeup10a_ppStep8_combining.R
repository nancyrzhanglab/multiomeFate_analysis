rm(list=ls())
library(Seurat)
library(Signac)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_data_rna.RData"))
all_data <- Seurat::CreateSeuratObject(counts = Seurat::GetAssayData(all_data_rna, layer = "counts"))

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
                           reduction.name = "Saver.pca")
print("Computing UMAP")
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            dims = 1:50,
                            reduction = "Saver.pca",
                            reduction.name = "Saver.umap")

#################

dataset_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM", "week5_CIS", "week5_COCL2", "week5_DABTRAM")
file_prefix <- "Writeup10a_ppStep7d_chromvar_"

for(dataset in dataset_vec){
  print(paste0("Working on ", dataset))
  load(paste0(out_folder, file_prefix, dataset, ".RData"))
  
  Seurat::DefaultAssay(all_data_subset) <- "chromvar"
  mat <- SeuratObject::LayerData(all_data_subset,
                                 layer = "data",
                                 assay = "chromvar")
  motif_list <- all_data_subset[["ATAC"]]@motifs@motif.names
  motif_names <- sapply(motif_names, function(x){x[1]})
  rownames(mat) <- motif_names
  
  mat_full <- Matrix::Matrix(0, 
                             nrow = nrow(mat), 
                             ncol = length(Seurat::Cells(all_data)), 
                             sparse = TRUE)
  rownames(mat_full) <- rownames(mat)
  colnames(mat_full) <- Seurat::Cells(all_data)
  mat_full[,colnames(mat)] <- mat
  
  all_data[[paste0("chromVar_", dataset)]] <- Seurat::CreateAssayObject(data = mat_full)
  all_data[[paste0("chromVar_", dataset)]]@var.features <- motif_list
}

#################

for(treatment in c("All", treatment_vec)){
  print(paste0("Working on ", treatment))
  
  mat <- read.csv(paste0(out_folder, "Writeup10a_ppStep7c_peakvi_", treatment, "_embedding.csv"))
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  mat <- as.matrix(mat)
  colnames(mat) <- paste0("peakVI", treatment, "_", 1:ncol(mat))
  
  set.seed(10)
  umap_res <- Seurat::RunUMAP(mat)
  umap_res <- umap_res@cell.embeddings
  
  mat_full <- Matrix::Matrix(0, 
                             nrow = length(Seurat::Cells(all_data)),
                             ncol = ncol(mat),
                             sparse = TRUE)
  rownames(mat_full) <- Seurat::Cells(all_data)
  colnames(mat_full) <- colnames(mat)
  mat_full[rownames(mat),] <- mat
  
  umap_full <- matrix(0, 
                      nrow = length(Seurat::Cells(all_data)),
                      ncol = 2)
  rownames(umap_full) <- Seurat::Cells(all_data)
  colnames(umap_full) <- colnames(umap_res)
  umap_full[rownames(umap_res),] <- umap_res
  
  all_data[[paste0("peakVI.", treatment)]] <- Seurat::CreateDimReducObject(mat)
  all_data[[paste0("pVI.", treatment, ".umap")]] <- Seurat::CreateDimReducObject(umap_full)
}

#################

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

all_data_saver <- all_data[["Saver"]]
all_data_saver_pca <- all_data[["Saver.pca"]]
all_data_saver_umap <- all_data[["Saver.umap"]]

all_data_fasttopic_CIS <- all_data[["fasttopic.CIS"]]
all_data_ft_CIS_umap <- all_data[["ft.CIS.umap"]]
all_data_fasttopic_COCL2 <- all_data[["fasttopic.COCL2"]]
all_data_ft_COCL2_umap <- all_data[["ft.COCL2.umap"]]
all_data_fasttopic_DABTRAM <- all_data[["fasttopic.DABTRAM"]]
all_data_ft_DABTRAM_umap <- all_data[["ft.DABTRAM.umap"]]
all_data_fasttopic_DABTRAM <- all_data[["fasttopic.DABTRAM"]]

all_data_peakVI_All <- all_data[["peakVI.All"]]
all_data_pVI_All_umap <- all_data[["pVI.All.umap"]]
all_data_peakVI_CIS <- all_data[["peakVI.CIS"]]
all_data_pVI_CIS_umap <- all_data[["pVI.CIS.umap"]]
all_data_peakVI_COCL2 <- all_data[["peakVI.COCL2"]]
all_data_pVI_COCL2_umap <- all_data[["pVI.COCL2.umap"]]
all_data_peakVI_DABTRAM <- all_data[["peakVI.DABTRAM"]]
all_data_pVI_DABTRAM_umap <- all_data[["pVI.DABTRAM.umap"]]

all_data_chromVar_day0 <- all_data[["chromVar_day0"]]
all_data_chromVar_day10_CIS <- all_data[["chromVar_day10_CIS"]]
all_data_chromVar_day10_COCL2 <- all_data[["chromVar_day10_COCL2"]]
all_data_chromVar_day10_DABTRAM <- all_data[["chromVar_day10_DABTRAM"]]
all_data_chromVar_week5_CIS <- all_data[["chromVar_week5_CIS"]]
all_data_chromVar_week5_COCL2 <- all_data[["chromVar_week5_COCL2"]]
all_data_chromVar_week5_DABTRAM <- all_data[["chromVar_week5_DABTRAM"]]

save(all_data_saver, all_data_saver_pca, all_data_saver_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver.RData"))
save(all_data_fasttopic_CIS, all_data_ft_CIS_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_fasttopic_CIS.RData"))
save(all_data_fasttopic_COCL2, all_data_ft_COCL2_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_fasttopic_COCL2.RData"))
save(all_data_fasttopic_DABTRAM, all_data_ft_DABTRAM_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_fasttopic_DABTRAM.RData"))

save(all_data_chromVar_day0, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_chromVar_day0.RData"))
save(all_data_chromVar_day10_CIS, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_chromVar_day10_CIS.RData"))
save(all_data_chromVar_day10_COCL2, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_chromVar_day10_COCL2.RData"))
save(all_data_chromVar_day10_DABTRAM, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_chromVar_day10_DABTRAM.RData"))
save(all_data_chromVar_week5_CIS, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_chromVar_week5_CIS.RData"))
save(all_data_chromVar_week5_COCL2, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_chromVar_week5_COCL2.RData"))
save(all_data_chromVar_week5_DABTRAM, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_chromVar_week5_DABTRAM.RData"))

save(all_data_peakVI_All, all_data_pVI_All_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_peakVI_All.RData"))
save(all_data_peakVI_CIS, all_data_pVI_CIS_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_peakVI_CIS.RData"))
save(all_data_peakVI_COCL2, all_data_pVI_COCL2_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_peakVI_COCL2.RData"))
save(all_data_peakVI_DABTRAM, all_data_pVI_DABTRAM_umap, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_peakVI_DABTRAM.RData"))

#################

# create an empty all_data object

mat <- Matrix::Matrix(0, 
                      ncol = length(Seurat::Cells(all_data)),
                      nrow = 1,
                      sparse = TRUE)
rownames(mat) <- "emptyFeature"
colnames(mat) <- Seurat::Cells(all_data)
all_data[["Empty"]] <- Seurat::CreateAssayObject(counts = mat)

Seurat::DefaultAssay(all_data) <- "Empty"
all_data <- Seurat::DietSeurat(all_data,
                               assays = "Empty")
load(paste0(out_folder, "Writeup10a_data_metadata.RData"))
var_names <- colnames(all_data_metadata)
var_names <- var_names[-grep("snn", var_names)]
var_names <- setdiff(var_names, "seurat_clusters")
all_data_metadata <- all_data_metadata[,var_names]
all_data@meta.data <- all_data_metadata

save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_empty.RData"))

print("Done! :)")

