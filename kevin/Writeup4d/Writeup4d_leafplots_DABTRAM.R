rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_DABTRAM.RData")

candidate_genes <- all_data_subset[["Saver"]]@var.features
candidate_genes <- intersect(candidate_genes, rownames(all_data_subset[["geneActivity"]]@data))
cor_vec <- sapply(1:length(candidate_genes), function(i){
  if(i %% floor(length(candidate_genes)/10) == 0) cat('*')
  
  gene <- candidate_genes[i]
  stats::cor(all_data_subset[["Saver"]]@data[gene,],
             all_data_subset[["geneActivity"]]@data[gene,])
})
cor_vec[which(is.na(cor_vec))] <- 0
quantile(cor_vec)
names(cor_vec) <- candidate_genes
cor_vec[order(cor_vec, decreasing = T)[1:18]]
selected_genes <- names(cor_vec)[order(cor_vec, decreasing = T)[1:18]]

datasets <- c("day0", "day10_DABTRAM", "week5_DABTRAM")
index_list <- lapply(datasets, function(x){
  idx <- which(all_data_subset$dataset == x)
  names(idx) <- NULL
  idx
})
names(index_list) <- datasets
col_vec <- rep(NA, ncol(all_data_subset))
for(i in 1:3){
  col_vec[index_list[[i]]] <- i
}

set.seed(10)
cell_idx <- sample(1:ncol(all_data_subset))
png("../../../../out/figures/Writeup4d/Writeup4d_DABTRAM_leafplots_RNA-ATAC.png", 
    height = 2500, width = 4500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(3,6))
for(gene in selected_genes){
  plot(y = all_data_subset[["geneActivity"]]@data[gene,cell_idx],
       x = all_data_subset[["Saver"]]@data[gene,cell_idx],
       ylab = "ATAC", xlab = "RNA", main = gene,
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
}
graphics.off()