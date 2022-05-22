rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_fasttopics_DABTRAM-geneactivity.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_DABTRAM.RData")

library(Seurat)
library(Signac)

all(rownames(topic_res$L) == colnames(all_data_subset))
all_data_subset[["fasttopic"]] <- Seurat::CreateDimReducObject(topic_res$L, 
                                                               assay = "RNA",
                                                               key = "fastTopic_")

geneactivity_mat <- all_data_subset[["geneActivity"]]@data
geneactivity_mat <- geneactivity_mat %*% all_data_subset@graphs[["Saver_snn"]] 

# Violin plots
Seurat::Idents(all_data_subset) <- "dataset"
plot1 <- Seurat::VlnPlot(all_data_subset, 
                         features = paste0("fastTopic_", 1:ncol(topic_res$L)),
                         ncol = 6,
                         pt.size = 0)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_geneactivity-DABTRAM_fasttopics-violin.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

round(apply(topic_res$F, 2, function(x){quantile(x, probs = c(0.75, 0.95, 0.99, 1))}),3)

idx <- order(topic_res$F[,17], decreasing = T)[1:18]
selected_genes <- rownames(topic_res$F)[idx]
selected_genes <- intersect(selected_genes, rownames(all_data_subset[["Saver"]]@data))
selected_genes <- intersect(selected_genes, rownames(all_data_subset[["unspliced"]]@data))

datasets <- c("day0", "day10_DABTRAM", "week5_DABTRAM")
index_list <- lapply(datasets, function(x){
  idx <- which(all_data_subset$dataset == x)
  names(idx) <- NULL
  idx
})
names(index_list) <- datasets
col_vec <- rep(NA, ncol(all_data_subset))
col_palette <- c(rgb(58, 120, 176, alpha = 255*0.25, maxColorValue = 255), 
                 rgb(240, 135, 55, alpha = 255*0.25, maxColorValue = 255), 
                 rgb(83, 157, 64, alpha = 255*0.25, maxColorValue = 255))
col_palette2 <- c(rgb(58, 120, 176, alpha = 255*0.8, maxColorValue = 255), 
                 rgb(240, 135, 55, alpha = 255*0.8, maxColorValue = 255), 
                 rgb(83, 157, 64, alpha = 255*0.8, maxColorValue = 255))
for(i in 1:3){
  col_vec[index_list[[i]]] <- col_palette[i]
}

set.seed(10)
cell_idx <- sample(1:ncol(all_data_subset))

png(paste0("../../../../out/figures/Writeup4d/Writeup4d_leafplots_DABTRAM_week5-ATACopen.png"), 
    height = 4500, width = 4500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(6,6))
for(i in 1:length(selected_genes)){
  cat('*')
  gene <- selected_genes[i]
  
  plot(y = all_data_subset[["unspliced"]]@data[gene,cell_idx],
       x = all_data_subset[["spliced"]]@data[gene,cell_idx],
       ylab = "Unspliced (Smoothed)", xlab = "Spliced (Smoothed)", 
       main = paste0(gene, ", Plot 1"),
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  
  x_means <- sapply(index_list, function(idx){
    mean(all_data_subset[["spliced"]]@data[gene,idx])
  })
  y_means <- sapply(index_list, function(idx){
    mean(all_data_subset[["unspliced"]]@data[gene,idx])
  })
  points(x = x_means, y = y_means, cex = 3.5, col = "white", pch = 16)
  points(x = x_means, y = y_means, cex = 2.5, col = col_palette2, pch = 16)
  
  axis(side = 1)
  axis(side = 2)
  
  plot(y = geneactivity_mat[gene,cell_idx],
       x = log1p(all_data_subset[["Saver"]]@data[gene,cell_idx]),
       ylab = "ATAC (Gene act., smoothed)", xlab = "RNA (SAVER, Log1p)", 
       main = paste0(gene, ", Plot 2"),
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  x_means <- sapply(index_list, function(idx){
    mean(log1p(all_data_subset[["Saver"]]@data[gene,idx]))
  })
  y_means <- sapply(index_list, function(idx){
    mean(geneactivity_mat[gene,idx])
  })
  points(x = x_means, y = y_means, cex = 3.5, col = "white", pch = 16)
  points(x = x_means, y = y_means, cex = 2.5, col = col_palette2, pch = 16)
  
  axis(side = 1)
  axis(side = 2)
}

graphics.off()

#############################

tmp <- topic_res$F
tmp <- tmp[-which(rownames(tmp) %in% selected_genes),]
tmp <- t(apply(tmp, 1, function(x){x/sum(x)}))
idx <- order(tmp[,17], decreasing = T)[1:18]
selected_genes2 <- rownames(tmp)[idx]

png(paste0("../../../../out/figures/Writeup4d/Writeup4d_leafplots_DABTRAM_week5-ATACopen2.png"), 
    height = 4500, width = 4500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(6,6))
for(i in 1:length(selected_genes2)){
  cat('*')
  gene <- selected_genes2[i]
  
  plot(y = all_data_subset[["unspliced"]]@data[gene,cell_idx],
       x = all_data_subset[["spliced"]]@data[gene,cell_idx],
       ylab = "Unspliced (Smoothed)", xlab = "Spliced (Smoothed)", 
       main = paste0(gene, ", Plot 1"),
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  
  x_means <- sapply(index_list, function(idx){
    mean(all_data_subset[["spliced"]]@data[gene,idx])
  })
  y_means <- sapply(index_list, function(idx){
    mean(all_data_subset[["unspliced"]]@data[gene,idx])
  })
  points(x = x_means, y = y_means, cex = 3.5, col = "white", pch = 16)
  points(x = x_means, y = y_means, cex = 2.5, col = col_palette2, pch = 16)
  
  axis(side = 1)
  axis(side = 2)
  
  plot(y = geneactivity_mat[gene,cell_idx],
       x = log1p(all_data_subset[["Saver"]]@data[gene,cell_idx]),
       ylab = "ATAC (Gene act., smoothed)", xlab = "RNA (SAVER, Log1p)", 
       main = paste0(gene, ", Plot 2"),
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  x_means <- sapply(index_list, function(idx){
    mean(log1p(all_data_subset[["Saver"]]@data[gene,idx]))
  })
  y_means <- sapply(index_list, function(idx){
    mean(geneactivity_mat[gene,idx])
  })
  points(x = x_means, y = y_means, cex = 3.5, col = "white", pch = 16)
  points(x = x_means, y = y_means, cex = 2.5, col = col_palette2, pch = 16)
  
  axis(side = 1)
  axis(side = 2)
}

graphics.off()

############################

tmp <- topic_res$F
tmp <- t(apply(tmp, 1, function(x){x/sum(x)}))
idx <- order(tmp[,19], decreasing = T)[1:18]
selected_genes3 <- rownames(tmp)[idx]
selected_genes3 <- intersect(selected_genes3, rownames(all_data_subset[["Saver"]]@data))
selected_genes3 <- intersect(selected_genes3, rownames(all_data_subset[["unspliced"]]@data))

png(paste0("../../../../out/figures/Writeup4d/Writeup4d_leafplots_DABTRAM_week5-ATACopen3.png"), 
    height = 4500, width = 4500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(6,6))
for(i in 1:length(selected_genes3)){
  cat('*')
  gene <- selected_genes3[i]
  
  plot(y = all_data_subset[["unspliced"]]@data[gene,cell_idx],
       x = all_data_subset[["spliced"]]@data[gene,cell_idx],
       ylab = "Unspliced (Smoothed)", xlab = "Spliced (Smoothed)", 
       main = paste0(gene, ", Plot 1"),
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  
  x_means <- sapply(index_list, function(idx){
    mean(all_data_subset[["spliced"]]@data[gene,idx])
  })
  y_means <- sapply(index_list, function(idx){
    mean(all_data_subset[["unspliced"]]@data[gene,idx])
  })
  points(x = x_means, y = y_means, cex = 3.5, col = "white", pch = 16)
  points(x = x_means, y = y_means, cex = 2.5, col = col_palette2, pch = 16)
  
  axis(side = 1)
  axis(side = 2)
  
  plot(y = geneactivity_mat[gene,cell_idx],
       x = log1p(all_data_subset[["Saver"]]@data[gene,cell_idx]),
       ylab = "ATAC (Gene act., smoothed)", xlab = "RNA (SAVER, Log1p)", 
       main = paste0(gene, ", Plot 2"),
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  x_means <- sapply(index_list, function(idx){
    mean(log1p(all_data_subset[["Saver"]]@data[gene,idx]))
  })
  y_means <- sapply(index_list, function(idx){
    mean(geneactivity_mat[gene,idx])
  })
  points(x = x_means, y = y_means, cex = 3.5, col = "white", pch = 16)
  points(x = x_means, y = y_means, cex = 2.5, col = col_palette2, pch = 16)
  
  axis(side = 1)
  axis(side = 2)
}

graphics.off()

