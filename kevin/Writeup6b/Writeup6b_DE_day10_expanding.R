rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
de_list <- vector("list", length = length(treatment_vec))
names(de_list) <- treatment_vec

Seurat::DefaultAssay(all_data) <- "RNA"
for(kk in 1:length(treatment_vec)){
  treatment <- treatment_vec[kk]
  print("====")
  print(treatment)
  
  ident_vec <- all_data$dataset
  names(ident_vec) <- colnames(all_data)
  lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
  cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
  cell_names2 <- colnames(all_data)[which(all_data$dataset == paste0("day10_", treatment))]
  cell_names_winning <- intersect(cell_names1, cell_names2)
  cell_names_losing <- setdiff(cell_names2, cell_names1)
  ident_vec[cell_names_winning] <- paste0("day10_win_", treatment)
  ident_vec[cell_names_losing] <- paste0("day10_lose_", treatment)
  all_data$ident <- ident_vec
  Seurat::Idents(all_data) <- "ident"
  
  set.seed(10)
  de_list[[kk]] <- Seurat::FindMarkers(all_data,
                                       ident.1 = paste0("day10_win_", treatment),
                                       ident.2 = paste0("day10_lose_", treatment),
                                       test.use = "wilcox",
                                       slot = "data",
                                       only.pos = F,
                                       verbose = F)
}

save(date_of_run, session_info, de_list,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_DE_day10_expanding.RData")


