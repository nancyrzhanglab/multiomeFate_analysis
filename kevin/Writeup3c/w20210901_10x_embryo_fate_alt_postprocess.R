rm(list=ls())
load("../../../../out/kevin/Writeup3c/w20210901_10x_mbrain_fate_alt_results.RData")

library(Seurat); library(Signac)

load_func <- function(){
  source("candidate_method_alt.R")
  source("chromatin_potential_alt.R")
  source("chromatin_potential_postprocess2.R")
  source("chromatin_potential_preparation_alt.R")
  source("diffusion_functions.R")
  source("estimation_method_alt.R")
  source("forming_method_alt.R")
  source("matching_method.R")
  source("metacell_construction.R")
  source("metacell_graph_plot.R")
}
load_func()

n <- nrow(diffusion_dist)
full_fate_prob <- matrix(0, n, ncol(fate_prob))
rownames(full_fate_prob) <- rownames(diffusion_dist)
for(i in 1:n){
  cell_name <- rownames(full_fate_prob)[i]
  if(cell_name %in% rownames(fate_prob)){
    idx <- which(rownames(fate_prob) == cell_name)
    full_fate_prob[i,] <- fate_prob[idx,]
  } else {
    stopifnot(cell_name %in% unlist(terminal_list))
    idx <- which(sapply(terminal_list, function(x){cell_name %in% x}))
    stopifnot(length(idx) == 1)
    full_fate_prob[i,idx] <- 1
  }
}
round(full_fate_prob, 2)

col_palette <- scales::hue_pal()(length(sort(unique(mbrain3@meta.data$celltype))))
col_gradient_list <- lapply(c("#FAE208", col_palette[5], "#5DFD35", col_palette[1]), function(col){
  grDevices::colorRampPalette(c("gray50", col))(10)
})
minmax_val <- min(apply(full_fate_prob,1,function(x){abs(diff(sort(x, decreasing = T)[1:2]))}))
seq_vec <- seq(minmax_val, 1, length.out = 10)
fate_col_vec <- sapply(1:nrow(full_fate_prob), function(i){
  idx <- which.max(full_fate_prob[i,])
  diff_val <- abs(diff(sort(full_fate_prob[i,], decreasing = T)[1:2]))
  bin <- which.min(abs(diff_val-seq_vec))
  col_gradient_list[[idx]][bin]
})

png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_fate_all.png"),
    height = 1300, width = 1400, res = 300, units = "px")
par(mar = c(3,3,4,0.5))
median_embedding <- compute_median_coords(mbrain3[["wnn.umap"]]@cell.embeddings, clustering)
plot_metacell_graph(median_embedding,
                    adj_mat,
                    col_vec = fate_col_vec,
                    xlim = c(min(median_embedding[,1]), 10),
                    xlab = "wnnUMAP_1",
                    ylab = "wnnUMAP_2",
                    main = "Fate probability to terminal cell-types",
                    axes = F,
                    xaxt = "n",
                    yaxt = "n",
                    asp = NA)
axis(side = 1)
axis(side = 2)
graphics.off()


# now plot the fates
target_name_vec <- c("Forebrain", "Oligo", "Cortical 2", "Cortical 1")
png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_fate.png"),
    height = 3000, width = 4500, res = 300, units = "px")
par(mfrow = c(2,3))
for(i in 1:ncol(fate_prob)){
  median_embedding <- compute_median_coords(mbrain3[["wnn.umap"]]@cell.embeddings, clustering)
  
  plot_metacell_graph(median_embedding,
                      adj_mat,
                      feature_vec = full_fate_prob[,i],
                      zlim = c(0.05,1),
                      xlab = "wnnUMAP_1",
                      ylab = "wnnUMAP_2",
                      main = paste0("Fate probability to ", target_name_vec[i]))
}
plot_legend(zlim = c(0.05,1),
            offset = -0.75)
graphics.off()

nonuniformity_vec <- apply(full_fate_prob, 1, function(x){
  r <- ncol(full_fate_prob)
  sum((x-1/r)^2)
})
nonuniformity_vec <- (nonuniformity_vec-min(nonuniformity_vec))/diff(range(nonuniformity_vec))

png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_nonuniformity.png"),
    height = 1300, width = 1400, res = 300, units = "px")
par(mar = c(3,3,4,0.5))
median_embedding <- compute_median_coords(mbrain3[["wnn.umap"]]@cell.embeddings, clustering)
plot_metacell_graph(median_embedding,
                    adj_mat,
                    feature_vec = nonuniformity_vec,
                    xlim = c(min(median_embedding[,1]), 10),
                    zlim = c(0,1),
                    xlab = "wnnUMAP_1",
                    ylab = "wnnUMAP_2",
                    main = "Fate deviation from uniformity",
                    axes = F,
                    xaxt = "n",
                    yaxt = "n",
                    asp = NA)
axis(side = 1)
axis(side = 2)
graphics.off()

#########################

gene_vec <- c("Pdgfra", "Sox10", "Cspg4", "Atp10b")
png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_oligo_genes.png"),
    height = 2000, width = 2000, res = 300, units = "px")
par(mfrow = c(2,2))

pred_mat <- .predict_yfromx2(mat_x, chrom_obj$res_g)
n <- nrow(mat_x)
col_vec <- rep("black", n)
col_vec[which(rownames(diffusion_dist) %in% terminal_list[[2]])] <- "red"
col_vec[which(rownames(diffusion_dist) %in% c("27", "45", "10", "15", "38"))] <- "green"
shuf <- unlist(lapply(c("black", "green", "red"), function(col){
  which(col_vec == col)
}))
col_vec2 <- col_vec[shuf]

for(gene in gene_vec){
  idx <- which(colnames(mat_y) == gene)
  x_vec <- mat_y[,idx]; y_vec <- pred_mat[,idx]
  x_vec <- x_vec[shuf]
  y_vec <- y_vec[shuf]
  plot(NA, 
       xlim = range(c(0,x_vec,y_vec)), 
       ylim = range(c(0,x_vec,y_vec)), 
       asp = T, 
       main = paste0("Oligo gene: ", gene),
       xlab = "RNA of meta-cell",
       ylab = "Predicted RNA from ATAC")
  lines(c(0,2*max(c(x_vec,y_vec))), 
        c(0,2*max(c(x_vec,y_vec))),
        col = "red",
        lwd = 2,
        lty = 2)
  points(x_vec,
         y_vec,
         col = col_vec2,
         pch = 16, cex = 2)
  
  if(gene == "Pdgfra"){
    legend("bottomright", 
           c("Oligo", "Glio", "Others"), 
           cex = 0.8, 
           fill = c("red", "green", "black"))
  }
}
graphics.off()

#########################

mat_x_full <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y_full <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))
peak_names <- sort(unique(mbrain3[["ATAC"]]@links$peak))
gene_names <- sort(unique(mbrain3[["ATAC"]]@links$gene))
mat_x_full <- mat_x_full[,which(colnames(mat_x_full) %in% peak_names)]
mat_y_full <- mat_y_full[,which(colnames(mat_y_full) %in% gene_names)]

gene_vec <- c("Pdgfra", "Sox10", "Cspg4", "Atp10b")
png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_oligo_genes_full.png"),
    height = 2000, width = 2000, res = 300, units = "px")
par(mfrow = c(2,2))

pred_mat <- .predict_yfromx2(mat_x_full, chrom_obj$res_g)
n <- nrow(mat_x_full)
col_vec <- rep("black", n)
col_vec[which(mbrain3@meta.data$celltype == "Oligodendrocyte")] <- "red"
col_vec[which(mbrain3@meta.data$celltype == "Glioblast")] <- "green"
shuf <- unlist(lapply(c("black", "green", "red"), function(col){
  which(col_vec == col)
}))
col_vec2 <- col_vec[shuf]

for(gene in gene_vec){
  idx <- which(colnames(mat_y_full) == gene)
  x_vec <- mat_y_full[,idx]; y_vec <- pred_mat[,idx]
  x_vec <- x_vec[shuf]
  y_vec <- y_vec[shuf]
  plot(NA, 
       xlim = range(c(0,x_vec,y_vec)), 
       ylim = range(c(0,x_vec,y_vec)), 
       asp = T, 
       main = paste0("Oligo gene: ", gene),
       xlab = "RNA of cells",
       ylab = "Predicted RNA from ATAC")
  lines(c(0,2*max(c(x_vec,y_vec))), 
        c(0,2*max(c(x_vec,y_vec))),
        col = "red",
        lwd = 2,
        lty = 2)
  points(x_vec,
         y_vec,
         col = col_vec2,
         pch = 16)
  
  if(gene == "Pdgfra"){
    legend("bottomright", 
           c("Oligo", "Glio", "Others"), 
           cex = 0.8, 
           fill = c("red", "green", "black"))
  }
}
graphics.off()

####################################


# first, plot all the meta-cells (no fate information)
png(paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_wnn_metacell_knn_colored.png"),
    height = 1500, width = 1400, res = 300, units = "px")
par(mar = c(3,3,4,0.5))
median_embedding <- compute_median_coords(mbrain3[["wnn.umap"]]@cell.embeddings, clustering)
col_palette <- scales::hue_pal()(length(sort(unique(mbrain3@meta.data$celltype))))
col_vec <- sapply(sort(unique(clustering)), function(x){
  celltype_vec <- mbrain3@meta.data$celltype[clustering == x]
  col_palette[which(sort(unique(mbrain3@meta.data$celltype)) == names(which.max(table(celltype_vec))))]
})
col_vec[which(rownames(adj_mat) %in% initial_vec)] <- col_palette[6]
plot_metacell_graph(median_embedding,
                    adj_mat,
                    col_vec = col_vec,
                    xlim = c(min(median_embedding[,1]), 10),
                    zlim = c(0,1),
                    xlab = "wnnUMAP_1",
                    ylab = "wnnUMAP_2",
                    main = "Mouse embryo (10x)\nSNN graph (Meta-cells by Seurat)",
                    axes = F,
                    xaxt = "n",
                    yaxt = "n",
                    asp = NA)
axis(side = 1)
axis(side = 2)
graphics.off()

# next, plot all the terminal cells
png(paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_wnn_metacell_terminal.png"),
    height = 1500, width = 1400, res = 300, units = "px")
par(mar = c(3,3,4,0.5))
median_embedding <- compute_median_coords(mbrain3[["wnn.umap"]]@cell.embeddings, clustering)
col_palette <- scales::hue_pal()(length(sort(unique(mbrain3@meta.data$celltype))))
col_vec <- sapply(sort(unique(clustering)), function(x){
  celltype_vec <- mbrain3@meta.data$celltype[clustering == x]
  main_celltype <- names(which.max(table(celltype_vec)))
  if(!x %in% c(unlist(terminal_list), initial_vec)) return("white")
  if(x %in% terminal_list[[3]]) return("#5DFD35")
  if(x %in% initial_vec) return(col_palette[6])
  col_palette[which(sort(unique(mbrain3@meta.data$celltype)) == main_celltype)]
})
plot_metacell_graph(median_embedding,
                    adj_mat,
                    col_vec = col_vec,
                    xlim = c(min(median_embedding[,1]), 10),
                    zlim = c(0,1),
                    xlab = "wnnUMAP_1",
                    ylab = "wnnUMAP_2",
                    main = "Initial/Terminal cell-types",
                    axes = F,
                    xaxt = "n",
                    yaxt = "n",
                    asp = NA)
axis(side = 1)
axis(side = 2)
graphics.off()

# now plot the cells based on the order they arrive into the method
for(iter in 1:max(chrom_obj$df_res$order_rec)){
  png(paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_umap_seurat_wnn_metacell_iter", iter, ".png"),
      height = 1500, width = 1400, res = 300, units = "px")
  par(mar = c(3,3,4,0.5))
  median_embedding <- compute_median_coords(mbrain3[["wnn.umap"]]@cell.embeddings, clustering)
  col_vec <- sapply(1:nrow(chrom_obj$df_res), function(i){
    if(chrom_obj$df_res$order_rec[i] > iter) return("white")
    if(chrom_obj$df_res$order_rec[i] == iter) return(grDevices::rgb(0.803, 0.156, 0.211))
    if(chrom_obj$df_res$order_rec[i] < iter) return("gray50")
  })
  plot_metacell_graph(median_embedding,
                      adj_mat,
                      col_vec = col_vec,
                      xlim = c(min(median_embedding[,1]), 10),
                      zlim = c(0,1),
                      xlab = "wnnUMAP_1",
                      ylab = "wnnUMAP_2",
                      main = paste0("Iteration ", iter),
                      axes = F,
                      xaxt = "n",
                      yaxt = "n",
                      asp = NA)
  axis(side = 1)
  axis(side = 2)
  graphics.off()
}
