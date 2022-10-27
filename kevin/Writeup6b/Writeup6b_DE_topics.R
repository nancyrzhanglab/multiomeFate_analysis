rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)

load("../../../../out/kevin/Writeup6/Writeup6_DE_day10-week5.RData")
load("../../../../out/kevin/Writeup6/Writeup6_all-data_lineage-assigned.RData")

de_week5_pairwise[["CIS"]][["combn_mat"]][,1:5]
head(de_week5_pairwise[["CIS"]][["de_list"]][[1]],20)
quantile(de_week5_pairwise[["CIS"]][["de_list"]][[1]][,"avg_log2FC"])

# first, find all lineages that really expand
treatment <- "week5_COCL2"
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
head(tab_mat[order(tab_mat[,treatment], decreasing = T),],10)
lineage_names <- rownames(tab_mat)[which(tab_mat[,treatment] >= 100)]

important_de_list <- sapply(lineage_names, function(lineage_name){
  idx <- which(de_week5_pairwise[["COCL2"]][["combn_mat"]] == lineage_name, arr.ind = T)[,2]
  tmp <- table(unlist(lapply(idx, function(kk){
    mat <- de_week5_pairwise[["COCL2"]][["de_list"]][[kk]]
    rownames(mat)[mat[,"p_val_adj"] <= 1e-4 & abs(mat[,"avg_log2FC"]) >= 1.5]
  })))
  sort(names(tmp)[which(tmp >= length(idx)/4)])
})
sapply(important_de_list, length)
names(important_de_list) <- lineage_names

################

for(i in 1:length(important_de_list)){
  print(paste0(i, " out of ", length(important_de_list)))
  if(length(important_de_list[[i]] > 1)){
    plot_topic_heatmap(topic_mat = all_data[["fasttopic_COCL2"]]@feature.loadings,
                       gene_vec = important_de_list[[i]],
                       file_name = paste0("../../../../out/figures/Writeup6b/Writeup6b_DE_COCL2_pairwiseWeek5_", names(important_de_list)[i], ".png"))
  }
}

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))] <- 1
table(keep_vec, all_data$dataset)
all_data$keep <- keep_vec
all_data_subset <- subset(all_data, keep == 1)

for(i in 1:length(important_de_list)){
  print(paste0(i, " out of ", length(important_de_list)))
  
  keep_vec <- all_data_subset$dataset
  keep_vec[all_data_subset$assigned_lineage == names(important_de_list)[i] & keep_vec == "week5_COCL2"] <- paste0("week5_", names(important_de_list)[i])
  all_data_subset$keep <- keep_vec
  
  Seurat::Idents(all_data_subset) <- "keep"
  plot1 <- Seurat::VlnPlot(all_data_subset, 
                           features = paste0("fastTopicCOCL2_", 1:ncol(all_data[["fasttopic_COCL2"]]@cell.embeddings)),
                           ncol = 5,
                           pt.size = 0)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6b/Writeup6b_COCL2_Week5_fasttopics-violin_", names(important_de_list)[i], ".png"),
                  plot1, device = "png", width = 15, height = 20, units = "in")
}
