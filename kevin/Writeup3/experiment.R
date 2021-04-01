rm(list=ls())
set.seed(10)
.rotate <- function(mat){t(mat)[,nrow(mat):1]}
p1 <- 100; p2 <- 20; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
set.seed(10)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window, 
                               signal_sd = 0.1)
timepoints <- 100; max_val <- 2
traj_mat <- generate_traj_cascading(df$df_y, timepoints = timepoints, 
                                    max_val = exp(max_val), min_val = 1)
set.seed(10)
obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list_traj_mat = list(traj_mat), verbose = T)
set.seed(10)
dat <- generate_data(obj_next, number_runs = 5, sample_perc = 1, time_tol = 0.01, 
                     verbose = T)
dim(dat$obs_x)

#############
vec_start <- which(dat$df_info$time <= 0.1)
list_end <- list(which(dat$df_info$time >= 0.9))
set.seed(10)
res <- chromatin_potential(dat$obs_x, dat$obs_y, df_x = dat$df_x, df_y = dat$df_y,
                           vec_start = vec_start, list_end = list_end, 
                           options = list(est_cis_window = window, cand_nn = 2))

image(.rotate(res$res_g$mat_g), asp = T)
quantile(res$res_g$mat_g)

set.seed(10)
time_vec <- dat$df_info$time + rnorm(n, sd = 0.01)
n <- nrow(dat$obs_x)
vec_from <- rep(NA, n); vec_to <- rep(NA, n)
key_vec <- as.character(sort(as.numeric(hash::keys(res$ht_neighbor))))
key_vec <- key_vec[order(res$df_res$order_rec, decreasing = F)]
for(i in 1:length(key_vec)){
  vec_from[i] <- time_vec[as.numeric(key_vec[i])]
  vec_to[i] <- time_vec[res$ht_neighbor[[key_vec[i]]][1]]
}
plot(NA, xlim = c(1,n), ylim = c(0,1))
for(i in 1:n){
  if(abs(vec_from[i] - vec_to[i]) <= 0.0001) {
    points(i, vec_from[i], pch = 16, col = 'red')
  } else {
    graphics::arrows(x0 = i, y0 = vec_from[i], x1 = i, y1 = vec_to[i], 
                     length = 0.1)
  }
}

####################

par(mfrow = c(1,2), mar = c(0.5, 0.5, 0.5, 0.5))
set.seed(10)
plot_umap(dat)
set.seed(10)
coords <- plot_umap(dat, only_coords = T)
plot(coords[,1], coords[,2], pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
key_vec <- hash::keys(res$ht_neighbor)
for(i in 1:n){
  idx1 <- as.numeric(key_vec[i])
  idx2 <- res$ht_neighbor[[key_vec[i]]][1]
  graphics::arrows(x0 = coords[idx1,1], y0 = coords[idx1,2],
                   x1 = coords[idx2,1], y1 = coords[idx2,2],
                   length = 0.05)
}
