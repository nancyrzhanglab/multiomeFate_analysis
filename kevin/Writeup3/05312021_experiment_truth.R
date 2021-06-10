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
                           x_sd_biological = 0.1, x_sd_technical = 0.5, 
                           y_exp_baseline = 0.1, y_sd_technical = 2,
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
dat <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)
idx <- which(df_cell$time >= 90)
dat$obs_x <- dat$obs_x[-idx,]; dat$obs_y <- dat$obs_y[-idx,]
dat$mean_x <- dat$mean_x[-idx,]; dat$mean_y <- dat$mean_y[-idx,]
df_cell <- df_cell[-idx,]

mat_x <- dat$mean_x; mat_y <- dat$mean_y

###############################

vec_start <- which(df_cell$time <= 10)
list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
  intersect(which(df_cell$branch == branch), which(df_cell$time >= 80))
})
# zz <- svd(mat_x); plot(zz$d[1:50])
rank_x <- 10
# zz <- svd(mat_y); plot(zz$d[1:50])
rank_y <- 10
set.seed(10)
prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                        vec_start, list_end,
                                        est_method = "threshold_glmnet",
                                        rec_method = "distant_cor_oracle",
                                        options = list(nn_nn = 20, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, est_cis_window = 30,
                                                       est_num_iterations = 4,
                                                       rec_bool_pred_nn = T))

set.seed(10)
res <- chromatin_potential(prep_obj, df_cell = df_cell, verbose = T, bool_oracle = T)
# set.seed(10)
# res2 <- chromatin_potential_refine(res, iter_max = 10)

save.image(file = "../../out/writeup3/05312021_true.RData")

#############

pred_y <- .predict_yfromx(mat_x, res$res_g, family = "gaussian")
png(file = "../../out/fig/writeup3/05312021_true_rna_predicted.png",
    height = 1500, width = 1500, res = 300, units = "px")
image(.rotate(pred_y))
graphics.off()
