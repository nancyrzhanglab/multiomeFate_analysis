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

res <- simulate_data_input(combn_wave_mat, 
                           x_exp_baseline = 0, x_exp_max = 0.7,
                           x_sd_biological = 0.5, x_sd_technical = 0.7, 
                           y_exp_baseline = 0.1, y_sd_technical = 2.5,
                           num_unrelated_x = 150, num_unrelated_y = 50, 
                           time_on = 15, time_windup = 15, 
                           max_lag = 25, min_lag = 10,
                           x_unrelated_intervals = 2,
                           x_unrelated_max = 1, y_unrelated_max = 4)
df_x <- res$df_x; df_y <- res$df_y
list_xnoise <- res$list_xnoise; list_ynoise <- res$list_ynoise

df_cell <- simulate_df_cell(5000, time_max = max(df_y$time_end_scaffold, na.rm = T),
                            num_branch = 3)

set.seed(10)
dat <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)
idx <- which(df_cell$time >= 90)
dat$obs_x <- dat$obs_x[-idx,]; dat$obs_y <- dat$obs_y[-idx,]
dat$mean_x <- dat$mean_x[-idx,]; dat$mean_y <- dat$mean_y[-idx,]
df_cell <- df_cell[-idx,]

mat_x <- dat$obs_x; mat_y <- dat$obs_y

###################

png(file = "../../out/fig/writeup3/06092021_data_atacrna_umap_colbybranch.png",
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,3))
set.seed(10)
zz <- Seurat::RunUMAP(mat_x)
idx <- sample(1:nrow(mat_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (ATAC)")

set.seed(10)
zz <- Seurat::RunUMAP(mat_y)
idx <- sample(1:nrow(mat_y))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (RNA)")

set.seed(10)
zz <- Seurat::RunUMAP(cbind(mat_x, mat_y))
idx <- sample(1:nrow(mat_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = (df_cell$branch+1)[idx], main = "UMAP (Both)")
graphics.off()

#####
col_palette <- colorRampPalette(c("red", "blue"))(10)
vec_val <- quantile(df_cell$time, probs = seq(0, 1, length.out = 10))
col_vec <- sapply(df_cell$time, function(x){col_palette[which.min(abs(vec_val - x))]})

png(file = "../../out/fig/writeup3/06092021_data_atacrna_umap_colbytime.png",
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,3))
set.seed(10)
zz <- Seurat::RunUMAP(mat_x)
idx <- sample(1:nrow(mat_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = col_vec[idx], main = "UMAP (ATAC)")

set.seed(10)
zz <- Seurat::RunUMAP(mat_y)
idx <- sample(1:nrow(mat_y))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = col_vec[idx], main = "UMAP (RNA)")

set.seed(10)
zz <- Seurat::RunUMAP(cbind(mat_x, mat_y))
idx <- sample(1:nrow(mat_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = col_vec[idx], main = "UMAP (Both)")
graphics.off()


png(file = "../../out/fig/writeup3/06092021_data_obs_atac.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(mat_x))
graphics.off()

png(file = "../../out/fig/writeup3/06092021_data_obs_rna.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(mat_y))
graphics.off()

