rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data) <- "chromvar"
treatmet_vec <- c("CIS", "COCL2", "DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

day10_lfc_list <- vector("list", length = 3)
names(day10_lfc_list) <- treatmet_vec
day10_pvalue_list <- day10_lfc_list
week5_lfc_list <- day10_lfc_list
week5_pvalue_list <- day10_lfc_list

# find the split in terms of Day10
for(treatment in treatmet_vec){
  print(paste0("Working on: Day10 of ", treatment))
  lineage_names_win <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
  cell_names_win <- colnames(all_data)[intersect(which(all_data$assigned_lineage %in% lineage_names_win),
                                                 which(all_data$dataset == paste0("day10_", treatment)))]
  lineage_names_lose <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] <= 5)]
  cell_names_lose <- colnames(all_data)[intersect(which(all_data$assigned_lineage %in% lineage_names_lose),
                                                  which(all_data$dataset == paste0("day10_", treatment)))]
  
  ident_vec <- rep(NA, ncol(all_data))
  names(ident_vec) <- colnames(all_data)
  ident_vec[cell_names_win] <- "winner"
  ident_vec[cell_names_lose] <- "loser"
  all_data$ident <- ident_vec
  Seurat::Idents(all_data) <- ident
  
  set.seed(10)
  de_res <- Seurat::FindMarkers(
    object = all_data,
    ident.1 = "winner",
    ident.2 = "loser",
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    verbose = F
  )
  print(head(de_res))
  
  day10_lfc_list[[treatment]] <- de_res$avg_diff
  day10_pvalue_list[[treatment]] <- de_res$p_val
}

# find the split in terms of Week5
for(treatment in treatmet_vec){
  print(paste0("Working on: Week5 of ", treatment))
  lineage_names_win <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 20)]
  cell_names_win <- colnames(all_data)[intersect(which(all_data$assigned_lineage %in% lineage_names_win),
                                                 which(all_data$dataset == paste0("week5_", treatment)))]
  lineage_names_lose <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 5)]
  cell_names_lose <- colnames(all_data)[intersect(which(all_data$assigned_lineage %in% lineage_names_lose),
                                                  which(all_data$dataset == paste0("week5_", treatment)))]
  
  ident_vec <- rep(NA, ncol(all_data))
  names(ident_vec) <- colnames(all_data)
  ident_vec[cell_names_win] <- "winner"
  ident_vec[cell_names_lose] <- "loser"
  all_data$ident <- ident_vec
  Seurat::Idents(all_data) <- ident
  
  set.seed(10)
  de_res <- Seurat::FindMarkers(
    object = all_data,
    ident.1 = "winner",
    ident.2 = "loser",
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    mean.fxn = rowMeans,
    fc.name = "avg_diff",
    verbose = F
  )
  print(head(de_res))
  
  week5_lfc_list[[treatment]] <- de_res$avg_diff
  week5_pvalue_list[[treatment]] <- de_res$p_val
}

save(date_of_run, session_info, 
     all_data, 
     day10_lfc_list, day10_pvalue_list,
     week5_lfc_list, week5_pvalue_list,
     file = "../../../../out/kevin/Writeup6n/Writeup6n_chromvar_de_all-timepoints.RData")
