rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_CIS.RData")

library(Seurat)
library(Signac)

# first find all the genes that qualify
gene_vec <- all_data_CIS[["RNA"]]@var.features
bool_vec <- sapply(gene_vec, function(gene){
  bool1 <- gene %in% rownames(all_data_CIS[["spliced"]]@data)
  bool2 <- gene %in% rownames(all_data_CIS[["maestro"]]@data)
  bool1 & bool2
})
candidate_genes <- sort(gene_vec[bool_vec])

########################

datasets <- c("day0", "day10_CIS", "week5_CIS")
index_list <- lapply(datasets, function(x){
  idx <- which(all_data_CIS$dataset == x)
  names(idx) <- NULL
  idx
})
names(index_list) <- datasets
  
# stress genes:
# ranking #1: are day-10 genes more spliced than day-0 genes
ranking1_val <- sapply(candidate_genes, function(gene){
  mean_day0 <- mean(all_data_CIS[["spliced"]]@data[gene,index_list[["day0"]]])
  sd_day0 <- sd(all_data_CIS[["spliced"]]@data[gene,index_list[["day0"]]])
  
  mean_day10 <- mean(all_data_CIS[["spliced"]]@data[gene,index_list[["day10_CIS"]]])
  sd_day10 <- sd(all_data_CIS[["spliced"]]@data[gene,index_list[["day10_CIS"]]])
  
  val <- (mean_day10 - mean_day0)/(sqrt(sd_day0^2+sd_day10^2))
  if(is.na(val)) val <- 0
  val
})
quantile(ranking1_val)
ranking1_vec <- rank(-ranking1_val)

# ranking #2: are week-5 genes less unspliced than day-10 genes
ranking2_val <- sapply(candidate_genes, function(gene){
  mean_day10 <- mean(all_data_CIS[["unspliced"]]@data[gene,index_list[["day10_CIS"]]])
  sd_day10 <- sd(all_data_CIS[["unspliced"]]@data[gene,index_list[["day10_CIS"]]])
  
  mean_week5 <- mean(all_data_CIS[["unspliced"]]@data[gene,index_list[["week5_CIS"]]])
  sd_week5 <- sd(all_data_CIS[["unspliced"]]@data[gene,index_list[["week5_CIS"]]])
  
  val <- (mean_day10 - mean_week5)/(sqrt(sd_week5^2+sd_day10^2))
  if(is.na(val)) val <- 0
  val
})
quantile(ranking2_val)
ranking2_vec <- rank(-ranking2_val)

# ranking #3: are day-10 genes more open than day-0 genes
ranking3_val <- sapply(candidate_genes, function(gene){
  mean_day0 <- mean(all_data_CIS[["maestro"]]@data[gene,index_list[["day0"]]])
  sd_day0 <- sd(all_data_CIS[["maestro"]]@data[gene,index_list[["day0"]]])
  
  mean_day10 <- mean(all_data_CIS[["maestro"]]@data[gene,index_list[["day10_CIS"]]])
  sd_day10 <- sd(all_data_CIS[["maestro"]]@data[gene,index_list[["day10_CIS"]]])
  
  val <- (mean_day10 - mean_day0)/(sqrt(sd_day0^2+sd_day10^2))
  if(is.na(val)) val <- 0
  val
})
quantile(ranking3_val)
ranking3_vec <- rank(-ranking3_val)

# ranking #4: are week-5 genes more open than day-0 genes
ranking4_val <- sapply(candidate_genes, function(gene){
  mean_day0 <- mean(all_data_CIS[["maestro"]]@data[gene,index_list[["day0"]]])
  sd_day0 <- sd(all_data_CIS[["maestro"]]@data[gene,index_list[["day0"]]])
  
  mean_week5 <- mean(all_data_CIS[["maestro"]]@data[gene,index_list[["week5_CIS"]]])
  sd_week5 <- sd(all_data_CIS[["maestro"]]@data[gene,index_list[["week5_CIS"]]])
  
  val <- (mean_week5 - mean_day0)/(sqrt(sd_week5^2+sd_day0^2))
  if(is.na(val)) val <- 0
  val
})
quantile(ranking4_val)
ranking4_vec <- rank(-ranking4_val)

######################

ranking_aggregate <- sapply(1:length(candidate_genes), function(i){
  max(c(ranking1_vec[i], ranking2_vec[i], ranking3_vec[i], ranking4_vec[i]))
})
names(ranking_aggregate) <- candidate_genes
idx <- order(ranking_aggregate, decreasing = F)[1:18]
selected_genes <- candidate_genes[idx]
col_vec <- rep(NA, ncol(all_data_CIS))
for(i in 1:3){
  col_vec[index_list[[i]]] <- i
}

png("../../../../out/figures/Writeup4d/Writeup4d_CIS_leafplots.png", 
    height = 4500, width = 4500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(6,6))
for(i in 1:length(selected_genes)){
  gene <- selected_genes[i]
  
  plot(y = all_data_CIS[["unspliced"]]@data[gene,],
       x = all_data_CIS[["spliced"]]@data[gene,],
       ylab = "Unspliced", xlab = "Spliced", main = selected_genes[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec, pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
  
  plot(y = all_data_CIS[["maestro"]]@data[gene,],
       x = all_data_CIS[["RNA"]]@scale.data[gene,],
       ylab = "ATAC", xlab = "RNA", main = selected_genes[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec, pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
}
graphics.off()

