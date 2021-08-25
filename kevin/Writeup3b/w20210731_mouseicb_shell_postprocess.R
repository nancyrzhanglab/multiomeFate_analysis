rm(list=ls())
load("../../../../out/kevin/Writeup3b/20210731_mouseicb_result.RData")

list_diagnos <- res$list_diagnos

.rowwise_median <- function(mat){
  vec <- matrixStats::colMedians(mat)
  dis <- sapply(1:nrow(mat), function(i){
    sum((mat[i,] - vec)^2)
  })
  mat[which.min(dis),]
}

png(file = "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_umap_timeoverlay.png",
    height = 3000, width = 3000, res = 300, units = "px")
mat_umap <- myeloid2[["umap.rna"]]@cell.embeddings
celltype <- as.factor(celltype)
col_palette <- scales::hue_pal()(length(levels(celltype)))
plot(mat_umap[,1], mat_umap[,2], asp = T, col = col_palette[as.numeric(celltype)], 
     pch = 16, main = "ATAC")
for(iter in 1:length(list_diagnos)){
  for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
    flip <- rbinom(1, 1, 0.1)
    if(flip == 1){
      idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
      vec_from <- .rowwise_median(mat_umap[idx_from,,drop=F])
      idx_to <- list_diagnos[[iter]]$recruit$postprocess[[i]]$to
      vec_to <- .rowwise_median(mat_umap[idx_to,,drop=F])
      
      graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                       x1 = vec_to[1], y1 = vec_to[2], length = 0.05)
    }
  }
}
graphics.off()

###################################

list_diagnos <- res$list_diagnos
n <- nrow(mat_x)
end_states <- c("Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")

# select endstate
prob_mat <- sapply(end_states, function(end_state){
  prob_vec <- rep(NA, n)
  prob_vec[which(celltype == end_state)] <- 1
  
  # set all other endstates to prob = 0
  prob_vec[which(celltype %in% end_states[which(end_states != end_state)])] <- 0
  
  # iterate backward
  while(any(is.na(prob_vec))){
    print(sum(is.na(prob_vec)))
    
    for(iter in 1:length(list_diagnos)){
      for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
        idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
        idx_to <- list_diagnos[[iter]]$recruit$postprocess[[i]]$to
        
        ## ignore all "another nodes" that don't have a value assigned yet
        valid_to <- idx_to[which(!is.na(prob_vec[idx_to]))]
        valid_from <- idx_from[which(is.na(prob_vec[idx_from]))]
        
        if(length(valid_from) == 0 | length(valid_to) == 0) next()
        ## do multiplication: prob of a node going into another node * prob of that node going to end state
        prob_vec[valid_from] <- mean(prob_vec[valid_to])
      }
    }
  }
  
  prob_vec
})



rownames(prob_mat) <- rownames(myeloid2@meta.data)
colnames(prob_mat) <- paste0("prob_", 1:ncol(prob_mat))
myeloid2[["prob"]] <- Seurat::CreateDimReducObject(embedding = prob_mat, key = "prob_")

plot1 <- Seurat::FeaturePlot(myeloid2, features = "prob_1", reduction = "wnn.umap")
plot1 <- plot1 + ggplot2::ggtitle("To Cluster 2")
plot2 <- Seurat::FeaturePlot(myeloid2, features = "prob_2", reduction = "wnn.umap")
plot2 <- plot2 + ggplot2::ggtitle("To Cluster 3")
plot3 <- Seurat::FeaturePlot(myeloid2, features = "prob_3", reduction = "wnn.umap")
plot3 <- plot3 + ggplot2::ggtitle("To Cluster 4")
plot4 <- Seurat::FeaturePlot(myeloid2, features = "prob_4", reduction = "wnn.umap")
plot4 <- plot4 + ggplot2::ggtitle("To Cluster 5")

tmp2 <- cowplot::plot_grid(plot1, plot2, plot3, plot4)
cowplot::save_plot(filename =  "../../../../out/figures/Writeup3b/Writeup3b_mouseicb_fateprob.png",
                   tmp2, ncol = 2, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")



