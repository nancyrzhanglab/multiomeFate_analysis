rm(list=ls())
g <- igraph::graph_from_edgelist(matrix(c(4,1, 4,5, 2,5, 3,5), nrow = 4, ncol = 2, byrow = T), 
                                 directed = F)
g <- igraph::set_vertex_attr(g, name = "lag", index = 4, value = 2)
g <- igraph::set_vertex_attr(g, name = "lag", index = 5, value = 5)
idx_root <- 4
num_waves <- 10
num_per_wave <- 5
distinct_waves <- 2

combn_wave_mat <- generate_combn_wave_mat(g, idx_root, num_waves = num_waves,
                                          num_per_wave = num_per_wave, 
                                          distinct_waves = distinct_waves)
res <- generate_data_input(combn_wave_mat, 
                           x_exp_baseline = 0.1, x_exp_max = 0.7,
                           x_sd_biological = 0.5, x_sd_technical = 0.2, 
                           y_exp_baseline = 0.1, y_sd_technical = 1,
                           num_unrelated_x = 100, num_unrelated_y = 50, 
                           time_on = 10, time_windup = 15, 
                           max_lag = 15, min_lag = 10,
                           x_unrelated_intervals = 2,
                           x_unrelated_max = 0.1, y_unrelated_max = 2)

# head(res$df_x)
# tail(res$df_x); tail(res$df_y)
# plot(res$df_x$time_start)
# plot(res$df_x$time_start, col = as.numeric(as.factor(res$df_x$branch)))
# plot(res$df_y$time_start_scaffold)
df_x <- res$df_x; df_y <- res$df_y
list_xnoise <- res$list_xnoise; list_ynoise <- res$list_ynoise

df_cell <- generate_df_cell(1000, time_max = max(df_x$time_end, na.rm = T),
                            num_branch = 3)

set.seed(10)
res <- generate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)

#####

image(.rotate(res$blueprint_x[[1]]))
image(.rotate(res$blueprint_x[[2]]))
image(.rotate(res$blueprint_x[[3]]))
image(.rotate(res$obs_x))
image(.rotate(res$obs_y))
image(.rotate(res$mean_x))
image(.rotate(res$mean_y))

col_palette <- colorRampPalette(c("red", "blue"))(10)
vec_val <- quantile(df_cell$time, probs = seq(0, 1, length.out = 10))
col_vec <- sapply(df_cell$time, function(x){col_palette[which.min(abs(vec_val - x))]})
gene_idx <- 80
peak_idx <- which(df_x$gene == df_y$name[gene_idx])
plot(rowSums(res$obs_x[,peak_idx]), res$obs_y[,gene_idx], pch = 16, col = col_vec)

gene_idx <- 80
peak_idx <- which(df_x$gene == df_y$name[gene_idx])
plot(df_cell$time, rowSums(res$obs_x[,peak_idx]), pch = 16, col = as.numeric(as.factor(df_cell$branch)))
plot(df_cell$time, rowSums(res$mean_x[,peak_idx]), pch = 16, col = as.numeric(as.factor(df_cell$branch)))
plot(df_cell$time, res$mean_y[,gene_idx], pch = 16, col = as.numeric(as.factor(df_cell$branch)))
plot(df_cell$time, res$obs_y[,gene_idx], pch = 16, col = as.numeric(as.factor(df_cell$branch)))
plot(df_cell$time, res$obs_y[,1], col = as.numeric(as.factor(df_cell$branch)), pch = 16)

par(mfrow = c(1,2))
set.seed(10)
zz <- Seurat::RunUMAP(res$obs_x)
idx <- sample(1:nrow(res$obs_x))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = as.numeric(as.factor(df_cell$branch))[idx],
     main = "UMAP (ATAC)")

set.seed(10)
zz <- Seurat::RunUMAP(res$obs_y)
idx <- sample(1:nrow(res$obs_y))
plot(zz@cell.embeddings[idx,1], zz@cell.embeddings[idx,2], asp = T, pch = 16,
     col = as.numeric(as.factor(df_cell$branch))[idx],
     main = "UMAP (RNA)")

################

col_palette <- colorRampPalette(c("red", "blue"))(10)
vec_val <- quantile(df_cell$time, probs = seq(0, 1, length.out = 10))
col_vec <- sapply(df_cell$time, function(x){col_palette[which.min(abs(vec_val - x))]})

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


     