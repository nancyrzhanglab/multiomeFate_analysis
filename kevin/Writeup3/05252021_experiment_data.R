rm(list=ls())
set.seed(10)
g <- igraph::graph_from_edgelist(matrix(c(4,1, 4,5, 2,5, 3,5), nrow = 4, ncol = 2, byrow = T), 
                                 directed = F)
g <- igraph::set_vertex_attr(g, name = "lag", index = 4, value = 3)
g <- igraph::set_vertex_attr(g, name = "lag", index = 5, value = 5)
idx_root <- 4; num_waves <- 10; num_per_wave <- 5; distinct_waves <- 2
combn_wave_mat <- simulate_combn_wave_mat(g, idx_root, num_waves = num_waves,
                                          num_per_wave = num_per_wave, 
                                          distinct_waves = distinct_waves)
plot(g)

res <- simulate_data_input(combn_wave_mat, 
                           x_exp_baseline = 0.1, x_exp_max = 0.7,
                           x_sd_biological = 0.1, x_sd_technical = 0.5, 
                           y_exp_baseline = 0.1, y_sd_technical = 2,
                           num_unrelated_x = 100, num_unrelated_y = 50, 
                           time_on = 15, time_windup = 15, 
                           max_lag = 20, min_lag = 10,
                           x_unrelated_intervals = 2,
                           x_unrelated_max = 0.1, y_unrelated_max = 2)
df_x <- res$df_x; df_y <- res$df_y
list_xnoise <- res$list_xnoise; list_ynoise <- res$list_ynoise
head(df_x); head(df_y)

df_cell <- simulate_df_cell(1000, time_max = max(df_y$time_end_scaffold, na.rm = T),
                            num_branch = 3)

set.seed(10)
res <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)
idx <- which(df_cell$time >= 90)
res$obs_x <- res$obs_x[-idx,]; res$obs_y <- res$obs_y[-idx,]
res$mean_x <- res$mean_x[-idx,]; res$mean_y <- res$mean_y[-idx,]
df_cell <- df_cell[-idx,]

#########################

png(file = "../../out/fig/writeup3/05252021_data_atac_blueprint.png",
    height = 1000, width = 2500, res = 300, units = "px")
par(mfrow = c(1,3))
for(i in 1:3){
  image(.rotate(res$blueprint_x[[i]]), main = paste0("Branch ", i), xlab = "Peak", ylab = "Cells")
  
}
graphics.off()

png(file = "../../out/fig/writeup3/05252021_data_atacrna_umap_colbybranch.png",
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,3))
set.seed(10)
zz <- Seurat::RunUMAP(res$obs_x )
idx <- sample(1:nrow(res$obs_x ))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (ATAC)")

set.seed(10)
zz <- Seurat::RunUMAP(res$obs_y)
idx <- sample(1:nrow(res$obs_y))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (RNA)")

set.seed(10)
zz <- Seurat::RunUMAP(cbind(res$obs_x, res$obs_y))
idx <- sample(1:nrow(res$obs_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (Both)")
graphics.off()

png(file = "../../out/fig/writeup3/05252021_data_atacrna_umap_true_colbybranch.png",
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,3))
set.seed(10)
zz <- Seurat::RunUMAP(res$mean_x)
idx <- sample(1:nrow(res$mean_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = as.numeric(as.factor(df_cell$branch))[idx], main = "UMAP (ATAC)")

set.seed(10)
zz <- Seurat::RunUMAP(res$mean_y)
idx <- sample(1:nrow(res$mean_y))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (RNA)")

set.seed(10)
zz <- Seurat::RunUMAP(cbind(res$mean_x, res$mean_y))
idx <- sample(1:nrow(res$mean_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (Both)")
graphics.off()

####

col_palette <- colorRampPalette(c("red", "blue"))(10)
vec_val <- quantile(df_cell$time, probs = seq(0, 1, length.out = 10))
col_vec <- sapply(df_cell$time, function(x){col_palette[which.min(abs(vec_val - x))]})

png(file = "../../out/fig/writeup3/05252021_data_atacrna_umap_colbytime.png",
    height = 1500, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2))
set.seed(10)
zz <- Seurat::RunUMAP(res$obs_x)
idx <- sample(1:nrow(res$obs_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = col_vec[idx], main = "UMAP (ATAC)")

set.seed(10)
zz <- Seurat::RunUMAP(res$obs_y)
idx <- sample(1:nrow(res$obs_y))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = col_vec[idx], main = "UMAP (RNA)")
graphics.off()

png(file = "../../out/fig/writeup3/05252021_data_atacrna_umap_true_colbytime.png",
    height = 1500, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2))
set.seed(10)
zz <- Seurat::RunUMAP(res$mean_x)
idx <- sample(1:nrow(res$mean_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = col_vec[idx], main = "UMAP (ATAC)")

set.seed(10)
zz <- Seurat::RunUMAP(res$mean_y)
idx <- sample(1:nrow(res$mean_y))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = col_vec[idx], main = "UMAP (RNA)")
graphics.off()

#####################################

png(file = "../../out/fig/writeup3/05252021_data_obs_atac.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(res$obs_x))
graphics.off()

png(file = "../../out/fig/writeup3/05252021_data_obs_rna.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(res$obs_y))
graphics.off()

png(file = "../../out/fig/writeup3/05252021_data_true_atac.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(res$mean_x))
graphics.off()

png(file = "../../out/fig/writeup3/05252021_data_true_rna.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(res$mean_y))
graphics.off()

#########################

col_palette <- colorRampPalette(c("red", "blue"))(10)
vec_val <- quantile(df_cell$time, probs = seq(0, 1, length.out = 10))
col_vec <- sapply(df_cell$time, function(x){col_palette[which.min(abs(vec_val - x))]})

for(gene_idx in c(1,49,111)){
  peak_idx <- which(df_x$gene == df_y$name[gene_idx])
  
  png(file = paste0("../../out/fig/writeup3/05252021_data_gene", gene_idx, ".png"),
      height = 2000, width = 3000, res = 300, units = "px")
  par(mfrow = c(2,4))
  plot(y = rowSums(res$obs_x[,peak_idx]), x = df_cell$time, col = df_cell$branch+1, pch = 16,
       xlab = "Time", ylab = "Expression", main = "ATAC (Observed)")
  plot(y = res$obs_y[,gene_idx], x = df_cell$time, col = df_cell$branch+1, pch = 16,
       xlab = "Time", ylab = "Expression", main = "RNA (Observed)")
  plot(x = rowSums(res$obs_x[,peak_idx]), y = res$obs_y[,gene_idx], col = df_cell$branch+1, 
       pch = 16,
       xlab = "ATAC (Obs)", ylab = "RNA (Obs)", main = "RNA vs. ATAC\n(By branch)")
  plot(x = rowSums(res$obs_x[,peak_idx]), y = res$obs_y[,gene_idx], col = col_vec, 
       pch = 16,
       xlab = "ATAC (Obs)", ylab = "RNA (Obs)", main = "RNA vs. ATAC\n(By time)")
  
  plot(y = rowSums(res$mean_x[,peak_idx]), x = df_cell$time, col = df_cell$branch+1, pch = 16,
       xlab = "Time", ylab = "Expression", main = "ATAC (True)")
  plot(y = res$mean_y[,gene_idx], x = df_cell$time, col = df_cell$branch+1, pch = 16,
       xlab = "Time", ylab = "Expression", main = "RNA (True)")
  plot(x = rowSums(res$mean_x[,peak_idx]), y = res$mean_y[,gene_idx], col = df_cell$branch+1, 
       pch = 16,
       xlab = "ATAC (True)", ylab = "RNA (True)", main = "RNA vs. ATAC\n(By branch)")
  plot(x = rowSums(res$mean_x[,peak_idx]), y = res$mean_y[,gene_idx], col = col_vec, 
       pch = 16,
       xlab = "ATAC (True)", ylab = "RNA (True)", main = "RNA vs. ATAC\n(By time)")
  graphics.off()
}

