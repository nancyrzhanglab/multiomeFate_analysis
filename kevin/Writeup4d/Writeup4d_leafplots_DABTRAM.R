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

# split into groups of 15
s_list <- lapply(1:ceiling(length(s.genes)/15), function(i){
  s.genes[((i-1)*15+1):(min(i*15, length(s.genes)))]
})
g2m_list <- lapply(1:ceiling(length(g2m.genes)/15), function(i){
  g2m.genes[((i-1)*15+1):(min(i*15, length(g2m.genes)))]
})

###############

geneactivity_mat <- all_data_subset[["geneActivity"]]@data
geneactivity_mat <- geneactivity_mat %*% all_data_subset@graphs[["Saver_snn"]] 

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
}

#############################################

candidate_genes <- all_data_subset[["Saver"]]@var.features
candidate_genes <-  intersect(candidate_genes, intersect(all_data_subset[["geneActivity"]]@var.features,
                                                         rownames(all_data_subset[["spliced"]])))

set.seed(10)
cell_subsample <- sample(1:ncol(all_data_subset), 5000)
cor_vec <- sapply(1:length(candidate_genes), function(i){
  print(i)
  
  gene <- candidate_genes[i]
  if(sum(all_data_subset[["unspliced"]]@data[gene,cell_subsample]) <= 1e-6 |
     sum(all_data_subset[["spliced"]]@data[gene,cell_subsample]) <= 1e-6) return(0)
  
  stats::cor(y = all_data_subset[["unspliced"]]@data[gene,cell_subsample],
             x = all_data_subset[["spliced"]]@data[gene,cell_subsample])
})
cor_vec[is.na(cor_vec)] <- 0
names(cor_vec) <- candidate_genes
selected_genes <- candidate_genes[order(cor_vec, decreasing = T)[1:18]]

png(paste0("../../../../out/figures/Writeup4d/Writeup4d_leafplots_DABTRAM_correlatedSplicedUnspliced.png"), 
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


set.seed(10)
cell_subsample <- sample(1:ncol(all_data_subset), 5000)
cor_vec2 <- sapply(1:length(candidate_genes), function(i){
  print(i)
  
  gene <- candidate_genes[i]
  if(sum(all_data_subset[["unspliced"]]@data[gene,cell_subsample]) <= 1e-6 |
     sum(all_data_subset[["spliced"]]@data[gene,cell_subsample]) <= 1e-6) return(0)
  
  stats::cor(y = geneactivity_mat[gene,cell_subsample],
             x = log1p(all_data_subset[["Saver"]]@data[gene,cell_subsample]))
})
cor_vec2[is.na(cor_vec2)] <- 0
names(cor_vec2) <- candidate_genes
selected_genes2 <- candidate_genes[order(cor_vec2, decreasing = T)[1:18]]

png(paste0("../../../../out/figures/Writeup4d/Writeup4d_leafplots_DABTRAM_correlatedATAC-RNA.png"), 
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
