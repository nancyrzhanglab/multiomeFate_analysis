rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_CIS.RData")
library(Seurat)
library(Signac)

Seurat::DefaultAssay(all_data_subset) <- "spliced"
all_data_subset[["spliced"]]@counts <- all_data_subset[["spliced"]]@data
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)

Seurat::DefaultAssay(all_data_subset) <- "unspliced"
all_data_subset[["unspliced"]]@counts <- all_data_subset[["unspliced"]]@data
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)

Seurat::DefaultAssay(all_data_subset) <- "maestro"
all_data_subset[["maestro"]]@counts <- all_data_subset[["maestro"]]@data
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)

##########################

tmp <- Matrix::rowSums(all_data_subset[["spliced"]]@counts)
remove_genes_spliced <- rownames(all_data_subset[["spliced"]])[which(tmp == 0)]
tmp <- Matrix::rowSums(all_data_subset[["unspliced"]]@counts)
remove_genes_unspliced <- rownames(all_data_subset[["unspliced"]])[which(tmp == 0)]
remove_genes <- unique(c(remove_genes_spliced, remove_genes_unspliced))

loading_mat <- all_data_subset[["fasttopic"]]@feature.loadings
quantile(log10(loading_mat[,"fastTopic_11"]))
vec <- loading_mat[,"fastTopic_11"]
vec <- vec[which(names(vec) %in% all_data_subset[["maestro"]]@var.features)]
vec <- vec[which(!names(vec) %in% remove_genes)]
quantile(log10(vec))
selected_genes <- names(vec)[order(vec, decreasing = T)[1:18]]

all(selected_genes %in% rownames(all_data_subset[["spliced"]]@scale.data))
all(selected_genes %in% rownames(all_data_subset[["unspliced"]]@scale.data))
all(selected_genes %in% rownames(all_data_subset[["maestro"]]@scale.data))
all(selected_genes %in% rownames(all_data_subset[["Saver"]]@data))

datasets <- c("day0", "day10_CIS", "week5_CIS")
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
png("../../../../out/figures/Writeup4d/Writeup4d_CIS_leafplots.png", 
    height = 4500, width = 4500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(6,6))
for(i in 1:length(selected_genes)){
  gene <- selected_genes[i]
  
  plot(y = all_data_subset[["unspliced"]]@scale.data[gene,cell_idx],
       x = all_data_subset[["spliced"]]@scale.data[gene,cell_idx],
       ylab = "Unspliced", xlab = "Spliced", main = selected_genes[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
  
  plot(y = all_data_subset[["maestro"]]@scale.data[gene,cell_idx],
       x = all_data_subset[["Saver"]]@data[gene,cell_idx],
       ylab = "ATAC", xlab = "RNA", main = selected_genes[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
}
graphics.off()

######################################################

load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_COCL2.RData")

Seurat::DefaultAssay(all_data_subset) <- "spliced"
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)

Seurat::DefaultAssay(all_data_subset) <- "unspliced"
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)

Seurat::DefaultAssay(all_data_subset) <- "maestro"
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)

tmp <- Matrix::rowSums(all_data_subset[["spliced"]]@counts)
remove_genes_spliced <- rownames(all_data_subset[["spliced"]])[which(tmp == 0)]
tmp <- Matrix::rowSums(all_data_subset[["unspliced"]]@counts)
remove_genes_unspliced <- rownames(all_data_subset[["unspliced"]])[which(tmp == 0)]
remove_genes <- unique(c(remove_genes_spliced, remove_genes_unspliced))

loading_mat <- all_data_subset[["fasttopic"]]@feature.loadings
quantile(log10(loading_mat[,"fastTopic_13"]))
vec <- loading_mat[,"fastTopic_13"]
vec <- vec[which(names(vec) %in% all_data_subset[["maestro"]]@var.features)]
vec <- vec[which(!names(vec) %in% remove_genes)]
quantile(log10(vec))
selected_genes <- names(vec)[order(vec, decreasing = T)[1:18]]

all(selected_genes %in% rownames(all_data_subset[["spliced"]]@scale.data))
all(selected_genes %in% rownames(all_data_subset[["unspliced"]]@scale.data))
all(selected_genes %in% rownames(all_data_subset[["maestro"]]@scale.data))
all(selected_genes %in% rownames(all_data_subset[["Saver"]]@data))

datasets <- c("day0", "day10_COCL2", "week5_COCL2")
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
png("../../../../out/figures/Writeup4d/Writeup4d_COCL2_leafplots.png", 
    height = 4500, width = 4500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(6,6))
for(i in 1:length(selected_genes)){
  gene <- selected_genes[i]
  
  plot(y = all_data_subset[["unspliced"]]@scale.data[gene,cell_idx],
       x = all_data_subset[["spliced"]]@scale.data[gene,cell_idx],
       ylab = "Unspliced", xlab = "Spliced", main = selected_genes[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
  
  plot(y = all_data_subset[["maestro"]]@scale.data[gene,cell_idx],
       x = all_data_subset[["Saver"]]@data[gene,cell_idx],
       ylab = "ATAC", xlab = "RNA", main = selected_genes[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
}
graphics.off()

selected_genes2 <- c("HSP90B1", "ABL2", "TRAF3", "HIST1H2AC", "IL7R", "PSMD13", "TNPO1", "SMYD4", "TIMM50")
selected_genes2 <- selected_genes2[which(selected_genes2 %in%  rownames(all_data_subset[["spliced"]]@scale.data))]
selected_genes2 <- selected_genes2[which(selected_genes2 %in%  rownames(all_data_subset[["maestro"]]@scale.data))]
selected_genes2 <- selected_genes2[which(selected_genes2 %in%  rownames(all_data_subset[["Saver"]]@scale.data))]

set.seed(10)
cell_idx <- sample(1:ncol(all_data_subset))
png("../../../../out/figures/Writeup4d/Writeup4d_COCL2_leafplots2.png", 
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = c(4, 4, 4, 0.5), mfrow = c(2,2))
for(i in 1:length(selected_genes2)){
  gene <- selected_genes2[i]
  
  plot(y = all_data_subset[["unspliced"]]@scale.data[gene,cell_idx],
       x = all_data_subset[["spliced"]]@scale.data[gene,cell_idx],
       ylab = "Unspliced", xlab = "Spliced", main = selected_genes2[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
  
  plot(y = all_data_subset[["maestro"]]@scale.data[gene,cell_idx],
       x = all_data_subset[["Saver"]]@data[gene,cell_idx],
       ylab = "ATAC", xlab = "RNA", main = selected_genes2[i],
       xaxt = "n", yaxt = "n", bty = "n",
       col = col_vec[cell_idx], pch = 16, cex = 1.5)
  axis(side = 1)
  axis(side = 2)
}
graphics.off()