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

##########################3

vec_start <- which(df_cell$time <= 10)
list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
  intersect(which(df_cell$branch == branch), which(df_cell$time >= 80))
})
# zz <- svd(mat_x); plot(zz$d[2:50])
rank_x <- 15
# zz <- svd(mat_y); plot(zz$d[2:50])
rank_y <- 15
set.seed(10)
prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                        vec_start, list_end,
                                        form_method = "average",
                                        est_method = "threshold_glmnet",
                                        rec_method = "distant_cor",
                                        options = list(nn_nn = 20, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, est_cis_window = 30,
                                                       est_num_iterations = 4,
                                                       rec_bool_pred_nn = T,
                                                       est_cv_choice = "lambda.min"))

set.seed(10)
res <- chromatin_potential(prep_obj,  verbose = T, bool_oracle = F)
save.image("../../out/writeup3/06092021_estimator.RData")

###################

n <- nrow(df_cell)
adj_mat <- matrix(0, n, n)
count_mat <- matrix(0, n, n)

for(i in 1:length(res$list_diagnos)){
  print(i)
  
  for(j in 1:length(res$list_diagnos[[i]]$recruit$postprocess)){
    from_idx <- res$list_diagnos[[i]]$recruit$postprocess[[j]]$from
    to_idx <- res$list_diagnos[[i]]$recruit$postprocess[[j]]$to
    
    avg_from <- as.numeric(.predict_yfromx(colMeans(mat_x[from_idx,,drop = F]), 
                                           res$res_g, family = "gaussian"))
    avg_to <- colMeans(mat_y[to_idx,,drop = F])
    val <- stats::cor(avg_from, avg_to, method = "kendall")
    
    adj_mat[from_idx, to_idx] <- adj_mat[from_idx, to_idx] + val
    count_mat[from_idx, to_idx] <- count_mat[from_idx, to_idx] + 1
  }
}

# apply a janky fix to the end-states: not all are reached
A <- adj_mat
S <- (A+t(A))/2
not_reached <- which(rowSums(S) == 0)
for(i in 1:3){
  to_idx <- intersect(not_reached, list_end[[i]])
  from_idx <- setdiff(list_end[[i]], not_reached)
  if(length(to_idx) == 0) next()
  
  avg_from <- as.numeric(.predict_yfromx(colMeans(mat_x[from_idx,,drop = F]), 
                                         res$res_g, family = "gaussian"))
  avg_to <- colMeans(mat_y[to_idx,,drop = F])
  val <- stats::cor(avg_from, avg_to, method = "kendall")
  
  adj_mat[from_idx, to_idx] <- adj_mat[from_idx, to_idx] + val
  count_mat[from_idx, to_idx] <- count_mat[from_idx, to_idx] + 1
}

count_mat <- pmax(count_mat, 1)
adj_mat <- adj_mat/count_mat
max_val <- max(adj_mat)

# now for the janky diffusion map
## see https://arxiv.org/pdf/1406.0013.pdf
A <- adj_mat
S <- (A+t(A))/2
q <- rowSums(S); q_inv <- 1/q
V <- .mult_mat_vec(.mult_vec_mat(q_inv, S), q_inv)
q1 <- rowSums(V); q1_inv <- 1/q1
Hss1 <- .mult_vec_mat(q1_inv, V)
eigen <- .svd_truncated(Hss1, 50, K_full_rank = F, vec_mean = NULL, vec_sd = NULL)
dm_coord <- eigen$u[,2:50]

# looks like garbage, let's abandon it

g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = T)
set.seed(10)
# coords <- igraph::layout_nicely(g)
# coords <- igraph::layout_with_lgl(g)
# coords <- igraph::layout_with_sugiyama(g) # took too long
# coords <- igraph::layout_with_mds(g)
coords <- ForceAtlas2::layout.forceatlas2(g, plotstep = 0)
plot(coords[,1], coords[,2], col = df_cell$branch+1, asp = T)

# all seem bad -- let's give up
set.seed(10)
obj <- Seurat::CreateSeuratObject(counts = Matrix::t(mat_x))
obj <- Seurat::ScaleData(obj)
obj <- Seurat::RunPCA(obj, features = rownames(obj), verbose = F)
obj <- Seurat::FindNeighbors(obj, dims = 1:15)
obj <- Seurat::FindClusters(obj, resolution = 0.5)
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

png(file = "../../out/fig/writeup3/06172021_est_atac_umap_timeoverlay.png",
    height = 3000, width = 3000, res = 300, units = "px")
plot(umap_res$embedding[,1], umap_res$embedding[,2], col = df_cell$branch+1, asp = T)

.l2norm <- function(x){sqrt(sum(x^2))}
for(i in 1:length(res$list_diagnos)){
  print(i)
  
  for(j in 1:length(res$list_diagnos[[i]]$recruit$postprocess)){
    samp <- stats::rbinom(1, 1, prob = 0.1)
    if(samp == 1){
      from_idx <- res$list_diagnos[[i]]$recruit$postprocess[[j]]$from
      to_idx <- res$list_diagnos[[i]]$recruit$postprocess[[j]]$to
      
      avg_from <- as.numeric(.predict_yfromx(colMeans(mat_x[from_idx,,drop = F]), 
                                             res$res_g, family = "gaussian"))
      avg_to <- colMeans(mat_y[to_idx,,drop = F])
      val <- stats::cor(avg_from, avg_to, method = "kendall")
      
      from_coord <- apply(umap_res$embedding[from_idx,], 2, median)
      to_coord <- apply(umap_res$embedding[to_idx,], 2, median)
      diff_vec <- to_coord - from_coord
      diff_vec <- val*diff_vec/.l2norm(diff_vec)/max_val
      
      graphics::arrows(x0 = from_coord[1], y0 = from_coord[2],
                       x1 = from_coord[1]+diff_vec[1], y1 = from_coord[2]+diff_vec[2], length = 0.05)
    }
  }
}
graphics.off()

# tmp <- assign_graph_attributes(g, df_cell)
# png(file = "../../out/fig/writeup3/06172021_est_nn.png",
#     height = 1500, width = 3000, res = 300, units = "px")
# par(mfrow = c(1,2), mar = rep(0.5, 4))
# set.seed(10)
# plot_igraph(tmp, color_by = "branch", bool_continuous = F, bool_default_par = F)
# set.seed(10)
# plot_igraph(tmp, color_by = "time", bool_continuous = T, bool_default_par = F)
# graphics.off()
