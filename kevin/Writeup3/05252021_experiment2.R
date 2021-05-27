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
                           x_exp_baseline = 0.1, x_exp_max = 0.7,
                           x_sd_biological = 0.5, x_sd_technical = 0.1, 
                           y_exp_baseline = 0.1, y_sd_technical = 1.5,
                           num_unrelated_x = 100, num_unrelated_y = 50, 
                           time_on = 15, time_windup = 15, 
                           max_lag = 20, min_lag = 10,
                           x_unrelated_intervals = 2,
                           x_unrelated_max = 0.1, y_unrelated_max = 2)
df_x <- res$df_x; df_y <- res$df_y
list_xnoise <- res$list_xnoise; list_ynoise <- res$list_ynoise

df_cell <- simulate_df_cell(1000, time_max = max(df_y$time_end_scaffold, na.rm = T),
                            num_branch = 3)

set.seed(10)
res <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell,
                     jitter = 0.01)
mat_x <- res$mean_x; mat_y <- res$mean_y

###############################

vec_start <- which(df_cell$time <= 10)
list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
  intersect(which(df_cell$branch == branch), which(df_cell$time >= 90))
})
# zz <- svd(mat_x); plot(zz$d[1:50])
rank_x <- 20
# zz <- svd(mat_y); plot(zz$d[1:50])
rank_y <- 20
set.seed(10)
prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                        vec_start, list_end,
                                        est_method = "threshold_glmnet",
                                        rec_method = "distant_cor_oracle",
                                        options = list(nn_nn = 10, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, est_cis_window = 30,
                                                       est_num_iterations = 4))

set.seed(10)
res <- chromatin_potential(prep_obj, df_cell = df_cell, verbose = T, bool_oracle = T)
set.seed(10)
res2 <- chromatin_potential_refine(res, iter_max = 10)

save.image(file = "../../out/writeup3/05252021_true.RData")

#############

png(file = "../../out/fig/writeup3/05252021_true_coefficient.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(res2$res_g$mat_g))
graphics.off()

plot(res2$res_g$vec_threshold, ylim = range(c(0,res2$res_g$vec_threshold)))
stopifnot(all(apply(pred_y, 2, min) >= res2$res_g$vec_threshold-1e-4))

pred_y <- .predict_yfromx(mat_x, res2$res_g, family = "gaussian")
png(file = "../../out/fig/writeup3/05252021_true_prediction_correlation.png",
    height = 1500, width = 1500, res = 300, units = "px")
xlim <- range(mat_y)
plot(NA, xlim = xlim, ylim = xlim, asp = T, xlab = "Predicted expression",
     ylab = "Obs expression")
for(i in 1:nrow(mat_x)){
  neigh <- res2$ht_neighbor[[as.character(i)]]
  obs_y <- colMeans(mat_y[neigh,,drop = F])
  points(x = pred_y[i,], y = obs_y, pch = 16, col = df_cell$branch[i])
}
graphics.off()

png(file = "../../out/fig/writeup3/05252021_true_prediction.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(pred_y), zlim = range(mat_y))
graphics.off()

png(file = "../../out/fig/writeup3/05252021_true_obs_x.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(mat_x))
graphics.off()

png(file = "../../out/fig/writeup3/05252021_true_obs_y.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(mat_y))
graphics.off()

##################

set.seed(10)
mat_umap <- Seurat::RunUMAP(mat_x)@cell.embeddings
plot(mat_umap[,1], mat_umap[,2], asp = T, col = df_cell$branch+1, pch = 16)
for(i in 1:nrow(mat_x)){
  flip <- rbinom(1, 1, 0.3)
  if(flip == 1){
    vec_from <- mat_umap[i,]
    idx_to <- res$ht_neighbor[[as.character(i)]]
    vec_to <- colMeans(mat_umap[idx_to,,drop=F])
    
    graphics::arrows(x0 = vec_from[1], y0 = vec_from[2],
                     x1 = vec_to[1], y1 = vec_to[2], length = 0.05)
  }
}

###

time_start <- df_cell$time
time_end <- sapply(1:nrow(mat_x), function(i){
  idx_to <- res$ht_neighbor[[as.character(i)]]
  mean(df_cell$time[idx_to])
})
time_start <- time_start[order(res$df_res$order_rec)]
time_end <- time_end[order(res$df_res$order_rec)]

plot(NA, xlim = c(0,nrow(mat_x)), ylim = range(time_start))
for(i in 1:length(time_start)){
  if(time_start[i] <= time_end[i]) col = "green" else col = "red"
  graphics::arrows(x0 = i, y0 = time_start[i],
                   x1 = i, y1 = time_end[i], length = 0.05, col = col)
}

