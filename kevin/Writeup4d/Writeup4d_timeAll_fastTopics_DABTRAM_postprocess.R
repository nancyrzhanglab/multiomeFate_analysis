rm(list=ls())

library(Seurat)
library(Signac)
library(fastTopics)

load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_peakmerging.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_fasttopics_DABTRAM.RData")

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_DABTRAM", "week5_DABTRAM"))] <- 1
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == 1)

topic_mat <- topic_res$L
topic_mat <- topic_mat[rownames(topic_mat) %in% colnames(all_data),]
topic_mat <- topic_mat[colnames(all_data),]

datasets <- c("day0", "day10_DABTRAM", "week5_DABTRAM")
index_list <- lapply(datasets, function(x){
  idx <- which(all_data$dataset == x)
  names(idx) <- NULL
  idx
})
names(index_list) <- datasets
col_vec <- rep(NA, ncol(all_data))
col_palette <- c("royalblue3", "orange2", "forestgreen")
for(i in 1:3){
  col_vec[index_list[[i]]] <- col_palette[i]
}

#####################################

set.seed(10)
cell_idx <- sample(1:ncol(all_data))
for(i in 1:ncol(topic_mat)){
  print(i)
  png(paste0("../../../../out/figures/Writeup4d/Writeup4d_topics_DABTRAM_scatterplot-",i, ".png"), 
      height = 3000, width = 3000, units = "px", res = 300)
  par(mfrow = c(5,6), mar = c(0.5, 0.5, 4, 0.5))
  for(j in 1:ncol(topic_mat)){
    plot(x = topic_mat[cell_idx,i],
         y = topic_mat[cell_idx,j],
         pch = 16, col = col_vec[cell_idx],
         main = paste0(i, " vs. ", j))
  }
  
  graphics.off()
}