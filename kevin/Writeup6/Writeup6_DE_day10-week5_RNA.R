rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6/Writeup6_all-data_lineage-assigned.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

head(table(all_data$assigned_lineage))
all(all_data$assigned_lineage[!is.na(all_data$assigned_lineage)] %in% rownames(all_data[["Lineage"]]@data))
quantile(all_data$assigned_posterior)
print("Finished loading in data")

###############################
## simply all Day10 vs Week5

Seurat::DefaultAssay(all_data) <- "Saver"
Seurat::Idents(all_data) <- all_data$dataset

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
set.seed(10)
de_day10vsweek5 <- lapply(treatment_vec, function(dataset){
  print(dataset)
  Seurat::FindMarkers(all_data,
                      ident.1 = paste0("day10_", dataset),
                      ident.2 = paste0("week5_", dataset),
                      test.use = "wilcox",
                      slot = "data",
                      only.pos = F,
                      verbose = F)
})
names(de_day10vsweek5) <- treatment_vec
print("Finished day10 vs week5")

###############################
## shrinking vs expanding, at day10

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
# tab_mat <- tab_mat[,c("day10_CIS", "week5_CIS")]
# tab_mat <- tab_mat[rowSums(tab_mat)>10,]
# tab_mat
# head(tab_mat[order(rowSums(tab_mat), decreasing = T),])

n <- ncol(all_data)

set.seed(10)
de_ExpandvsShrunk <- lapply(treatment_vec, function(dataset){
  print(dataset)
  tab_mat2 <- tab_mat[,c(paste0("day10_", dataset), paste0("week5_", dataset))]
  tab_mat2 <- tab_mat2[rowSums(tab_mat2)>=10,]
  print(head(tab_mat2))
  
  key_vec <- rep(NA, n)
  for(i in 1:nrow(tab_mat2)){
    if(tab_mat2[i,paste0("week5_", dataset)] >= 10) {
      key_vec[which(all_data$assigned_lineage == rownames(tab_mat2)[i] & all_data$dataset == paste0("day10_", dataset))] <- "expand"
    } else {
      key_vec[which(all_data$assigned_lineage == rownames(tab_mat2)[i] & all_data$dataset == paste0("day10_", dataset))] <- "shrunk"
    }
  }
  all_data$key <- key_vec
  Seurat::Idents(all_data) <- all_data$key
  
  Seurat::FindMarkers(all_data,
                      ident.1 = "expand",
                      ident.2 = "shrunk",
                      test.use = "wilcox",
                      slot = "data",
                      only.pos = F,
                      verbose = F)
})
names(de_ExpandvsShrunk) <- treatment_vec
print("Finished expand vs shrink, day10")

# head(de_ExpandvsShrunk)
# de_ExpandvsShrunk[de_ExpandvsShrunk$p_val <= 1e-3,]
# dim(de_ExpandvsShrunk)
# quantile(de_ExpandvsShrunk$p_val)

###############################
## among the different expanding

de_week5_pairwise <- lapply(treatment_vec, function(dataset){
  print(dataset)
  lineage_expanding <- rownames(tab_mat)[tab_mat[,paste0("week5_", dataset)] >= 10]
  
  n <- ncol(all_data)
  key_vec <- rep(NA, n)
  for(i in 1:length(lineage_expanding)){
    lineage_name <- lineage_expanding[i]
    key_vec[which(all_data$assigned_lineage == lineage_name & all_data$dataset == paste0("day10_", dataset))] <- paste0(lineage_name, "_day10")
    key_vec[which(all_data$assigned_lineage == lineage_name & all_data$dataset == paste0("week5_", dataset))] <- paste0(lineage_name, "_week5")
  }
  table(key_vec)
  all_data$key <- key_vec
  Seurat::Idents(all_data) <- all_data$key
  
  combn_mat <- utils::combn(length(lineage_expanding), 2)
  tmp1 <- lineage_expanding[combn_mat[1,]]
  tmp2 <- lineage_expanding[combn_mat[2,]]
  combn_mat <- rbind(tmp1, tmp2)
  
  de_list <- lapply(1:ncol(combn_mat), function(j){
    print(paste0(j, " out of ", ncol(combn_mat)))
    
    set.seed(10)
    Seurat::FindMarkers(all_data,
                        ident.1 = paste0(combn_mat[1,j], "_week5"),
                        ident.2 = paste0(combn_mat[2,j], "_week5"),
                        test.use = "wilcox",
                        slot = "data",
                        only.pos = F,
                        verbose = F)
  })
  
  list(de_list = de_list, combn_mat = combn_mat)
})
names(de_week5_pairwise) <- treatment_vec
print("Finished week5 expand")


save(date_of_run, session_info,
     de_day10vsweek5, de_ExpandvsShrunk,
     de_week5_pairwise, tab_mat,
     file = "../../../../out/kevin/Writeup6/Writeup6_DE_day10-week5.RData")

# sapply(de_week5_pairwise, function(x){
#   length(which(x$p_val <= 1e-4))
# })
# vec <- table(unlist(lapply(de_week5_pairwise, function(x){
#   rownames(x)[x$p_val <= 1e-4]
# })))
# vec[order(vec)]
# quantile(vec)

de_day10_pairwise <- lapply(treatment_vec, function(dataset){
  print(dataset)
  lineage_expanding <- rownames(tab_mat)[tab_mat[,paste0("week5_", dataset)] >= 10]
  
  n <- ncol(all_data)
  key_vec <- rep(NA, n)
  for(i in 1:length(lineage_expanding)){
    lineage_name <- lineage_expanding[i]
    key_vec[which(all_data$assigned_lineage == lineage_name & all_data$dataset == paste0("day10_", dataset))] <- paste0(lineage_name, "_day10")
    key_vec[which(all_data$assigned_lineage == lineage_name & all_data$dataset == paste0("week5_", dataset))] <- paste0(lineage_name, "_week5")
  }
  table(key_vec)
  all_data$key <- key_vec
  Seurat::Idents(all_data) <- all_data$key
  
  combn_mat <- utils::combn(length(lineage_expanding), 2)
  tmp1 <- lineage_expanding[combn_mat[1,]]
  tmp2 <- lineage_expanding[combn_mat[2,]]
  combn_mat <- rbind(tmp1, tmp2)
  
  de_list <- lapply(1:ncol(combn_mat), function(j){
    print(paste0(j, " out of ", ncol(combn_mat)))
    
    ident_1 <- paste0(combn_mat[1,j], "_day10")
    ident_2 <- paste0(combn_mat[2,j], "_day10")
    if(length(which(all_data$key == ident_1)) <= 3 || length(which(all_data$key == ident_2)) <= 3){
      tmp <- matrix(NA, nrow = 0, ncol = 5)
      colnames(tmp) <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
      return(tmp)
    } else {
      set.seed(10)
      Seurat::FindMarkers(all_data,
                          ident.1 = ident_1,
                          ident.2 = ident_2,
                          test.use = "wilcox",
                          slot = "data",
                          only.pos = F,
                          verbose = F)
    }
  })
  
  list(de_list = de_list, combn_mat = combn_mat)
})
names(de_day10_pairwise) <- treatment_vec
print("Finished day10 expand")

# quantile(sapply(de_day10_pairwise, function(x){
#   length(which(x$p_val <= 1e-4))
# }))
# vec <- table(unlist(lapply(de_day10_pairwise, function(x){
#   rownames(x)[x$p_val <= 1e-4]
# })))
# vec[order(vec)]
# quantile(vec)


save(date_of_run, session_info,
     de_day10vsweek5, de_ExpandvsShrunk,
     de_week5_pairwise, de_day10_pairwise,
     tab_mat,
     file = "../../../../out/kevin/Writeup6/Writeup6_DE_day10-week5.RData")



