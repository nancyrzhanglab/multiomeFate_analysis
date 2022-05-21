rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_DABTRAM.RData")

library(Seurat); library(Signac)

s.genes <- cc.genes$s.genes
s.genes <- intersect(s.genes, intersect(all_data_subset[["Saver"]]@var.features, 
                                        intersect(all_data_subset[["geneActivity"]]@var.features,
                                                  rownames(all_data_subset[["spliced"]]))))
s.genes <- sort(s.genes)
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- intersect(g2m.genes, intersect(all_data_subset[["Saver"]]@var.features, 
                                        intersect(all_data_subset[["geneActivity"]]@var.features,
                                                  rownames(all_data_subset[["spliced"]]))))
g2m.genes <- sort(g2m.genes)

datasets <- c("day0", "day10_DABTRAM", "week5_DABTRAM")
index_list <- lapply(datasets, function(x){
  idx <- which(all_data_subset$dataset == x)
  names(idx) <- NULL
  idx
})
names(index_list) <- datasets
col_vec <- rep(NA, ncol(all_data_subset))
col_palette <- c(rgb(58, 120, 176, maxColorValue = 255), 
                 rgb(240, 135, 55, maxColorValue = 255), 
                 rgb(83, 157, 64, maxColorValue = 255))
for(i in 1:3){
  col_vec[index_list[[i]]] <- col_palette[i]
}

set.seed(10)
cell_idx <- sample(1:ncol(all_data_subset))

# split into groups of 15
s_list <- lapply(1:ceiling(length(s.genes)/15), function(i){
  s.genes[((i-1)*15+1):(min(i*15, length(s.genes)))]
})
g2m_list <- lapply(1:ceiling(length(g2m.genes)/15), function(i){
  g2m.genes[((i-1)*15+1):(min(i*15, length(g2m.genes)))]
})

###############


for(panel_idx in 1:length(s_list)){
  panel <- s_list[[panel_idx]]
  print(panel_idx)
  
  png(paste0("../../../../out/figures/Writeup4d/Writeup4d_leafplots_DABTRAM_cellcycle_s-", panel_idx, ".png"), 
      height = 4500, width = 4500, units = "px", res = 300)
  par(mar = c(4, 4, 4, 0.5), mfrow = c(5,6))
  for(i in 1:length(panel)){
    cat('*')
    gene <- panel[i]
    
    plot(y = all_data_subset[["unspliced"]]@data[gene,cell_idx],
         x = all_data_subset[["spliced"]]@data[gene,cell_idx],
         ylab = "Unspliced", xlab = "Spliced", main = gene,
         xaxt = "n", yaxt = "n", bty = "n",
         col = col_vec[cell_idx], pch = 16, cex = 1.5)
    axis(side = 1)
    axis(side = 2)
    
    plot(y = all_data_subset[["geneActivity"]]@data[gene,cell_idx],
         x = all_data_subset[["Saver"]]@data[gene,cell_idx],
         ylab = "ATAC", xlab = "RNA", main = gene,
         xaxt = "n", yaxt = "n", bty = "n",
         col = col_vec[cell_idx], pch = 16, cex = 1.5)
    axis(side = 1)
    axis(side = 2)
  }
  
  graphics.off()
}

for(panel_idx in 1:length(g2m_list)){
  panel <- g2m_list[[panel_idx]]
  print(panel_idx)
  
  png(paste0("../../../../out/figures/Writeup4d/Writeup4d_leafplots_DABTRAM_cellcycle_g2m-", panel_idx, ".png"), 
      height = 4500, width = 4500, units = "px", res = 300)
  par(mar = c(4, 4, 4, 0.5), mfrow = c(5,6))
  for(i in 1:length(panel)){
    cat('*')
    gene <- panel[i]
    
    plot(y = all_data_subset[["unspliced"]]@data[gene,cell_idx],
         x = all_data_subset[["spliced"]]@data[gene,cell_idx],
         ylab = "Unspliced", xlab = "Spliced", main = gene,
         xaxt = "n", yaxt = "n", bty = "n",
         col = col_vec[cell_idx], pch = 16, cex = 1.5)
    axis(side = 1)
    axis(side = 2)
    
    plot(y = all_data_subset[["geneActivity"]]@data[gene,cell_idx],
         x = all_data_subset[["Saver"]]@data[gene,cell_idx],
         ylab = "ATAC", xlab = "RNA", main = gene,
         xaxt = "n", yaxt = "n", bty = "n",
         col = col_vec[cell_idx], pch = 16, cex = 1.5)
    axis(side = 1)
    axis(side = 2)
  }
  
  graphics.off()
}
