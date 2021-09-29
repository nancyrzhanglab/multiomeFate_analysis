rm(list=ls())
load("../../../../out/kevin/Writeup3c/w20210901_10x_mbrain_fate_alt_results.RData")

df_res <- chrom_obj$df_res
P <- chrom_obj$snn
n <- nrow(P)

for(i in 1:n){
  idx <- which(P[i,] != 0)
  min_val <- min(P[i,idx])
  max_val <- max(P[i,idx])
  P[i,idx] <- exp(-(P[i,idx]- min_val)/max_val)
}
P <- diag(1/matrixStats::rowSums2(P)) %*% P

r <- max(df_res$init_state, na.rm = T)
terminal_list <- lapply(1:r, function(k){
  which(df_res$init_state == k)
})
terminal_vec <- unlist(terminal_list)

A <- solve(diag(n - length(terminal_vec)) - P[-terminal_vec, -terminal_vec]) %*% P[-terminal_vec, terminal_vec]
colnames(A) <- terminal_vec
rownames(A) <- colnames(diffusion_dist)[-terminal_vec]

fate_prob <- sapply(1:r, function(k){
  idx <- which(colnames(A) %in% as.character(terminal_list[[k]]))
  rowSums(A[,idx,drop = F])
})
colnames(fate_prob) <- as.character(1:r)

##################

initial_vec <- c("4", "48")
terminal_list <- list(c("28", "36", "19", "11"), #forebrain
                      c("46"), # oligo
                      c("52", "6", "0", "8"), #cortical2
                      c("1", "2", "3", "7", #cortical1
                        "14", "18", 
                        "24", "25", "29", 
                        "30", "32", "37", "39",
                        "40", "42","44", "47",
                        "53"))
diffusion_dist <- chrom_obj$diffusion_dist

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

png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_fate_nosignal_all.png"),
    height = 1300, width = 1400, res = 300, units = "px")
par(mar = c(3,3,4,0.5))
median_embedding <- compute_median_coords(mbrain3[["wnn.umap"]]@cell.embeddings, clustering)
plot_metacell_graph(median_embedding,
                    adj_mat,
                    col_vec = fate_col_vec,
                    xlim = c(min(median_embedding[,1]), 10),
                    xlab = "wnnUMAP_1",
                    ylab = "wnnUMAP_2",
                    main = "Fate probability to terminal cell-types\n(Only based on one modality)",
                    axes = F,
                    xaxt = "n",
                    yaxt = "n",
                    asp = NA)
axis(side = 1)
axis(side = 2)
graphics.off()


target_name_vec <- c("Forebrain", "Oligo", "Cortical 2", "Cortical 1")
png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_fate_nosignal.png"),
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
                      main = paste0("Fate probability to ", target_name_vec[i],"\nNo direction"))
}
plot_legend(zlim = c(0,1),
            offset = -0.75)
graphics.off()

nonuniformity_vec <- apply(full_fate_prob, 1, function(x){
  r <- ncol(full_fate_prob)
  sum((x-1/r)^2)
})
nonuniformity_vec <- (nonuniformity_vec-min(nonuniformity_vec))/diff(range(nonuniformity_vec))

png(paste0("../../../../out/figures/Writeup3c/w20210901_10x_embryo_nonuniformity_nosignal.png"),
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
                    main = "Fate deviation from uniformity\n(Only based on one modality)",
                    axes = F,
                    xaxt = "n",
                    yaxt = "n",
                    asp = NA)
axis(side = 1)
axis(side = 2)
graphics.off()

#################

col_palette <- scales::hue_pal()(length(sort(unique(mbrain3@meta.data$celltype))))
col_gradient_list <- lapply(c(col_palette[2], "#FAE208", "#5DFD35", col_palette[1]), function(col){
  grDevices::colorRampPalette(c("gray50", col))(10)
})
for(i in 1:length(col_gradient_list)){
  png(paste0("../../../../out/figures/Writeup3c/legend_", i, ".png"),
      height = 1500, width = 800, res = 300, units = "px")
  plot_legend(zlim = c(0,1),
              high_color = col_gradient_list[[i]][10],
              low_color = col_gradient_list[[i]][1])
  graphics.off()
}
