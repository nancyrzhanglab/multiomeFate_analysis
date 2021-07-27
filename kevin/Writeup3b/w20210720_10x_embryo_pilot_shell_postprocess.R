rm(list=ls())
load("../../../../out/kevin/Writeup3b/10x_embryo.RData")
load("../../../../out/kevin/Writeup3b/10x_embryo_result_tmp.RData")

# list_diagnos <- res$list_diagnos

.rowwise_median <- function(mat){
  vec <- matrixStats::colMedians(mat)
  dis <- sapply(1:nrow(mat), function(i){
    sum((mat[i,] - vec)^2)
  })
  mat[which.min(dis),]
}

png(file = "../../../../out/figures/Writeup3b/Writeup3b_10x_embryo_umap_timeoverlay2.png",
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
for(iter in 1:length(list_diagnos)){
  for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
    flip <- rbinom(1, 1, 0.5)
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

#######
set.seed(10)
tmp <- mat_y; tmp <- scale(tmp, center = T, scale = F)
pc_res <- irlba::irlba(tmp, nv = 50)
tmp2 <- multiomeFate:::.mult_mat_vec(pc_res$u, pc_res$d)
mat_umap <- Seurat::RunUMAP(tmp2)@cell.embeddings
celltype <- as.factor(celltype)
col_palette <- scales::hue_pal()(length(levels(celltype)))
for(iter in 1:length(list_diagnos)){
  png(file = paste0("../../../../out/figures/Writeup3b/Writeup3b_10x_embryo_umap_timeoverlay2_iter", iter, ".png"),
      height = 3000, width = 3000, res = 300, units = "px")
  plot(mat_umap[,1], mat_umap[,2], asp = T, col = col_palette[as.numeric(celltype)], 
       pch = 16, main = "ATAC")
  for(i in 1:length(list_diagnos[[iter]]$recruit$postprocess)){
    flip <- rbinom(1, 1, 0.9)
    if(flip == 1){
      idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
      vec_from <- .rowwise_median(mat_umap[idx_from,,drop=F])
      idx_to <- list_diagnos[[iter]]$recruit$postprocess[[i]]$to
      vec_to <- .rowwise_median(mat_umap[idx_to,,drop=F])
      
      graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                       x1 = vec_to[1], y1 = vec_to[2], length = 0.05)
    }
  }
  graphics.off()
}

####################

iter <- 10
zz <- sapply(1:length(list_diagnos[[iter]]$recruit$postprocess), function(i){
  idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
  print(i)
  print(table(celltype[idx_from]))
  matrixStats::colMedians(mat_umap[idx_from,,drop=F])
})
t(zz)
tmp_idx <- 3
idx_to <- list_diagnos[[iter]]$recruit$postprocess[[tmp_idx]]$to
matrixStats::colMedians(mat_umap[idx_to,,drop=F])

i <- 27
idx_from <- list_diagnos[[iter]]$recruit$postprocess[[i]]$from
data.frame(idx_from, celltype[idx_from])

####################

mbrain2 <- Seurat::CreateSeuratObject(counts = t(mat_x))
mbrain2[["celltype"]] <- celltype
rownames(mat_umap) <- rownames(mat_x)
mbrain2[["asdf"]] <- Seurat::CreateDimReducObject(embedding = mat_umap, assay = "RNA")

plot1 <- Seurat::DimPlot(mbrain2, reduction = "asdf", 
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("10x Mouse Embryo\nRNA") 
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3b/Writeup3b_10x_embryo_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


