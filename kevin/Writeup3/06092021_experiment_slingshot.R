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

##################

set.seed(10)
obj <- Seurat::CreateSeuratObject(counts = Matrix::t(mat_x))
obj <- Seurat::ScaleData(obj)
obj <- Seurat::RunPCA(obj, features = rownames(obj), verbose = F)
obj <- Seurat::FindNeighbors(obj, dims = 1:15)
obj <- Seurat::FindClusters(obj, resolution = 0.5)
set.seed(10)
obj <- Seurat::RunUMAP(obj, dims = 1:15)

plot1 <- Seurat::DimPlot(obj, reduction = "umap")
plot1 <- plot1 + ggplot2::ggtitle("ATAC (Seurat default clustering)")
ggplot2::ggsave(filename = "../../out/fig/writeup3/06092021_slingshot_atac_clustering.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

set.seed(10)
tmp <- Seurat::RunUMAP(mat_x)@cell.embeddings
obj[["umap2"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP2")
plot1 <- Seurat::DimPlot(obj, reduction = "umap2")
plot1 <- plot1 + ggplot2::ggtitle("ATAC (Seurat default clustering)")
ggplot2::ggsave(filename = "../../out/fig/writeup3/06092021_slingshot_atac_clustering2.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###########

cluster_vec <- as.numeric(obj@meta.data$seurat_clusters)
#modify clusterings by the true starts and ends
vec_start <- which(df_cell$time <= 10)
list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
  intersect(which(df_cell$branch == branch), which(df_cell$time >= 80))
})
num_clust <- max(cluster_vec)
cluster_vec[vec_start] <- num_clust+1
for(i in 1:length(list_end)){
  cluster_vec[list_end[[i]]] <- num_clust+1+i
}
table(cluster_vec)


start_clus_str <- as.character(num_clust+1)
end_clus_str <- sapply(num_clust+1+1:length(list_end), function(x){as.character(x)})
dimmat_x <- obj[["pca"]]@cell.embeddings[,1:15]
set.seed(10)
slingshot_lin <- slingshot::getLineages(dimmat_x, cluster_vec, start.clus = start_clus_str, 
                    end.clus = end_clus_str)
slingshot_curve <- slingshot::getCurves(slingshot_lin)

## add the slingshot points onto the plot
set.seed(10)
umap_res <- uwot::umap(obj[["pca"]]@cell.embeddings[,1:15],
                       metric = "cosine", n_epochs = NULL,
                       learning_rate = 1.0,
                       min_dist = 0.3,
                       spread = 1.0,
                       set_op_mix_ratio = 1.0,
                       local_connectivity = 1L,
                       repulsion_strength = 1,
                       negative_sample_rate = 5,
                       a = NULL,
                       b = NULL,
                       fast_sgd = F,
                       ret_nn = T, ret_model = T, verbose = T)
curve_points <- do.call(rbind, lapply(1:3, function(x){
  slingshot_curve@curves[[x]]$s[slingshot_curve@curves[[x]]$ord,]
}))
set.seed(10)
umap_res2 <- uwot::umap_transform(X = curve_points, model = umap_res)

png(file = paste0("../../out/fig/writeup3/06092021_slingshot_trajectory.png"),
    height = 3000, width = 3000, res = 300, units = "px")
plot(umap_res$embedding[,1], umap_res$embedding[,2], pch = 16,
     cex = 2)
graphics::points(umap_res$embedding[,1], umap_res$embedding[,2], pch = 16,
                 col = df_cell$branch+1)
cutoff <- nrow(umap_res2)/3
for(i in 1:3){
  idx <-(((i-1)*cutoff+1):(i*cutoff))
  graphics::lines(umap_res2[idx,1], umap_res2[idx,2],
                  lwd = 5, col = "white")
  graphics::lines(umap_res2[idx,1], umap_res2[idx,2],
                  lwd = 1, col = i+1)
}
graphics.off()


