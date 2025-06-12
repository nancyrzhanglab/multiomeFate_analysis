rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"

all_data <- multiomeFate::data_loader(which_files = c("saver"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  keep_vec <- rep(FALSE, ncol(all_data))
  keep_vec[which(all_data$dataset %in% c("day0", 
                                         paste0("day10_", treatment), 
                                         paste0("week5_", treatment)))] <- TRUE
  all_data$keep <- keep_vec
  all_data2 <- subset(all_data, keep == TRUE)
  
  Seurat::DefaultAssay(all_data2) <- "Saver"
  all_data2 <- Seurat::ScaleData(all_data2)
  all_data2 <- Seurat::RunPCA(all_data2, 
                             verbose = FALSE,
                             reduction.name = paste0("Saver.", treatment, ".pca"),
                             reduction.key = paste0("Saver", treatment, "PC_"))
  
  mat <- all_data2[[paste0("Saver.", treatment, ".pca")]]@cell.embeddings
  cell_names <- Seurat::Cells(all_data)
  mat_full <- matrix(NA, nrow = length(cell_names), ncol = ncol(mat))
  colnames(mat_full) <- colnames(mat)
  rownames(mat_full) <- cell_names
  mat_full[rownames(mat),] <- mat
  all_data[[paste0("Saver.", treatment, ".pca")]] <- Seurat::CreateDimReducObject(mat_full)
}

dataset_vec <- unique(all_data$dataset)

for(dataset in dataset_vec){
  print(paste0("Working on ", dataset))
  all_data2 <- subset(all_data, dataset == dataset)
  
  Seurat::DefaultAssay(all_data2) <- "Saver"
  all_data2 <- Seurat::NormalizeData(all_data2)
  all_data2 <- Seurat::FindVariableFeatures(all_data2)
  all_data2 <- Seurat::ScaleData(all_data2)
  all_data2 <- Seurat::RunPCA(all_data2, 
                              verbose = FALSE,
                              reduction.name = paste0("Saver.", dataset, ".pca"),
                              reduction.key = paste0("Saver", dataset, "PC_"))
  
  mat <- all_data2[[paste0("Saver.", dataset, ".pca")]]@cell.embeddings
  cell_names <- Seurat::Cells(all_data)
  mat_full <- matrix(NA, nrow = length(cell_names), ncol = ncol(mat))
  colnames(mat_full) <- colnames(mat)
  rownames(mat_full) <- cell_names
  mat_full[rownames(mat),] <- mat
  all_data[[paste0("Saver.", dataset, ".pca")]] <- Seurat::CreateDimReducObject(mat_full)
}

all_data_saver_CIS_pca <- all_data[["Saver.CIS.pca"]]
all_data_saver_COCL2_pca <- all_data[["Saver.COCL2.pca"]]
all_data_saver_DABTRAM_pca <- all_data[["Saver.DABTRAM.pca"]]

all_data_saver_day0_pca <- all_data[["Saver.day0.pca"]]
all_data_saver_day10_CIS_pca <- all_data[["Saver.day10_CIS.pca"]]
all_data_saver_day10_COCL2_pca <- all_data[["Saver.day10_COCL2.pca"]]
all_data_saver_day10_DABTRAM_pca <- all_data[["Saver.day10_DABTRAM.pca"]]
all_data_saver_week5_CIS_pca <- all_data[["Saver.week5_CIS.pca"]]
all_data_saver_week5_COCL2_pca <- all_data[["Saver.week5_COCL2.pca"]]
all_data_saver_week5_DABTRAM_pca <- all_data[["Saver.week5_DABTRAM.pca"]]

save(all_data_saver_CIS_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_CIS_pca.RData"))
save(all_data_saver_COCL2_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_COCL2_pca.RData"))
save(all_data_saver_DABTRAM_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_DABTRAM_pca.RData"))

save(all_data_saver_day0_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_day0_pca.RData"))
save(all_data_saver_day10_CIS_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_day10_CIS_pca.RData"))
save(all_data_saver_day10_COCL2_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_day10_COCL2_pca.RData"))
save(all_data_saver_day10_DABTRAM_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_day10_DABTRAM_pca.RData"))
save(all_data_saver_week5_CIS_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_week5_CIS_pca.RData"))
save(all_data_saver_week5_COCL2_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_week5_COCL2_pca.RData"))
save(all_data_saver_week5_DABTRAM_pca, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_saver_week5_DABTRAM_pca.RData"))

print("Done! :)")
