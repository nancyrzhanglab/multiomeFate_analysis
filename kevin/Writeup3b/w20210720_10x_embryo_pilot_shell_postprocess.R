rm(list=ls())
load("../../../../out/kevin/Writeup3b/10x_embryo_result.RData")

png(file = "../../../../out/figures/Writeup3b/Writeup3b_10x_embryo_umap_timeoverlay.png",
    height = 3000, width = 3000, res = 300, units = "px")
set.seed(10)
tmp <- mat_y; tmp <- scale(tmp, center = T, scale = F)
pc_res <- irlba::irlba(tmp, nv = 50)
tmp2 <- multiomeFate:::.mult_mat_vec(pc_res$u, pc_res$d)
mat_umap <- Seurat::RunUMAP(tmp2)@cell.embeddings
celltype <- as.factor(celltype)
col_palette <- scales::hue_pal()(length(levels(celltype)))
plot(mat_umap[,1], mat_umap[,2], asp = T, col = col_palette[as.numeric(celltype)], 
     pch = 16, main = "ATAC")
for(iter in 1:length(res$list_diagnos)){
  for(i in 1:length(res$list_diagnos[[iter]]$recruit$postprocess)){
    flip <- rbinom(1, 1, 0.1)
    if(flip == 1){
      idx_from <- res$list_diagnos[[iter]]$recruit$postprocess[[i]]$from
      vec_from <- matrixStats::colMedians(mat_umap[idx_from,,drop=F])
      idx_to <- res$list_diagnos[[iter]]$recruit$postprocess[[i]]$to
      vec_to <- matrixStats::colMedians(mat_umap[idx_to,,drop=F])
      
      graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                       x1 = vec_to[1], y1 = vec_to[2], length = 0.05)
    }
  }
}
graphics.off()

#######
set.seed(10)
tmp <- mat_y; tmp <- scale(tmp, center = T, scale = F)
pc_res <- irlba::irlba(tmp, nv = 50)
tmp2 <- multiomeFate:::.mult_mat_vec(pc_res$u, pc_res$d)
mat_umap <- Seurat::RunUMAP(tmp2)@cell.embeddings
celltype <- as.factor(celltype)
col_palette <- scales::hue_pal()(length(levels(celltype)))
for(iter in 1:length(res$list_diagnos)){
  png(file = paste0("../../../../out/figures/Writeup3b/Writeup3b_10x_embryo_umap_timeoverlay_iter", iter, ".png"),
      height = 3000, width = 3000, res = 300, units = "px")
  plot(mat_umap[,1], mat_umap[,2], asp = T, col = col_palette[as.numeric(celltype)], 
       pch = 16, main = "ATAC")
  for(i in 1:length(res$list_diagnos[[iter]]$recruit$postprocess)){
    flip <- rbinom(1, 1, 0.9)
    if(flip == 1){
      idx_from <- res$list_diagnos[[iter]]$recruit$postprocess[[i]]$from
      vec_from <- matrixStats::colMedians(mat_umap[idx_from,,drop=F])
      idx_to <- res$list_diagnos[[iter]]$recruit$postprocess[[i]]$to
      vec_to <- matrixStats::colMedians(mat_umap[idx_to,,drop=F])
      
      graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                       x1 = vec_to[1], y1 = vec_to[2], length = 0.05)
    }
  }
  graphics.off()
}

