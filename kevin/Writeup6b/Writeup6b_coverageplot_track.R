rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("gene_list.R")
source("coverage_extractor.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::Idents(all_data) <- "dataset"
Seurat::DefaultAssay(all_data) <- "ATAC"

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
gene_vec <- sort(unique(unlist(keygenes)))
gene_vec <- gene_vec[which(gene_vec %in% rownames(all_data[["RNA"]]))]

ident_mat <- sapply(treatment_vec, function(treatment){
  ident_vec <- all_data$dataset
  names(ident_vec) <- colnames(all_data)
  lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
  cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
  cell_names2 <- colnames(all_data)[which(all_data$dataset == "day0")]
  cell_names_winning <- intersect(cell_names1, cell_names2)
  cell_names_losing <- setdiff(cell_names2, cell_names1)
  ident_vec[cell_names_winning] <- paste0("day0_win_", treatment)
  ident_vec[cell_names_losing] <- "day0_lose"
  ident_vec
})
colnames(ident_mat) <- paste0("day0_win_", treatment_vec)

# finalize the idents, mainly to find the day0 cells that lost in all 3 trials
ident_vec <- apply(ident_mat, 1, function(x){
  uniq_vec <- unique(x)
  if(length(uniq_vec) == 1) {
    return(uniq_vec[1])
  } else return("N/A")
})
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"

track_list_1 <- lapply(gene_vec, function(gene){
  print(gene)
  coverage_extractor(
    object = all_data,
    region = gene,
    group.by = "ident", 
    which_ident = c("day0_lose", "day10_CIS", "day10_COCL2", "day10_DABTRAM", "week5_CIS", "week5_COCL2", "week5_DABTRAM"),
    assay = "ATAC",
    extend.downstream = 1000,
    extend.upstream = 1000,
    window = 100)
})

save(track_list_1, session_info, date_of_run,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_coverageplot_track-tmp.RData")

print("Moving onto the second phase")
# now compute all the other tracklist for the day0 winners
track_ll <- vector("list", 3)
names(track_ll) <- colnames(ident_mat)
for(kk in 1:3){
  print("=============")
  print(kk)
  print("=============")
  all_data$ident <- ident_mat[,kk]
  Seurat::Idents(all_data) <- "ident"
  
  track_ll[[kk]] <- lapply(gene_vec, function(gene){
    print(gene)
    coverage_extractor(
      object = all_data,
      region = gene,
      group.by = "ident", 
      which_ident = colnames(ident_mat)[kk],
      assay = "ATAC",
      extend.downstream = 1000,
      extend.upstream = 1000,
      window = 100)
  })
}

save(track_list_1, track_ll, session_info, date_of_run,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_coverageplot_track-tmp.RData")

###################

# do some merging

track_list <- track_list_1
len <- length(track_list)
for(j in 1:len){
  for(kk in 1:3){
    track_list[[j]]$coverage_mean <- cbind(track_list[[j]]$coverage_mean, track_ll[[kk]][[j]]$coverage_mean[,1])
    colnames(track_list[[j]]$coverage_mean)[ncol(track_list[[j]]$coverage_mean)] <- colnames(track_ll[[kk]][[j]]$coverage_mean)[1]
  }
  
  track_list[[j]]$coverage_mean <- track_list[[j]]$coverage_mean[,colnames(track_list[[j]]$coverage_mean)]
 
  for(kk in 1:3){
    track_list[[j]]$coverage_count <- cbind(track_list[[j]]$coverage_count, track_ll[[kk]][[j]]$coverage_count)
    colnames(track_list[[j]]$coverage_count)[ncol(track_list[[j]]$coverage_count)] <- colnames(track_ll[[kk]][[j]]$coverage_mean)[1]
  }
  
  track_list[[j]]$coverage_count <- track_list[[j]]$coverage_count[,colnames(track_list[[j]]$coverage_count)]
  
  for(kk in 1:3){
    track_list[[j]]$total_vec <- c(track_list[[j]]$total_vec, track_ll[[kk]][[j]]$total_vec[1])
    names(track_list[[j]]$total_vec)[length(track_list[[j]]$total_vec)] <- colnames(track_ll[[kk]][[j]]$coverage_mean)[1]
  }
  
  track_list[[j]]$total_vec <- track_list[[j]]$total_vec[names(track_list[[j]]$total_vec)]
}

names(track_list) <- gene_vec

save(track_list, session_info, date_of_run,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_coverageplot_track.RData")


