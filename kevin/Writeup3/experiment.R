rm(list=ls())
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
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
dat <- generate_data(obj_next, number_runs = 10, sample_perc = 1, time_tol = 0.01, 
                     verbose = T)
dim(dat$obs_x)

##################

vec_start <- which(dat$df_info$time <= 0.1) # these are the cells at the start state
list_end <- list(which(dat$df_info$time >= 0.9)) # these are the cells at the end state
set.seed(10)
# run the estimator
res <- chromatin_potential(dat$obs_x, dat$obs_y, df_x = dat$df_x, df_y = dat$df_y,
                           vec_start = vec_start, list_end = list_end, 
                           mat_g_init = mat_g, 
                           form_method = "literal", est_method = "glmnet",
                           cand_method = "nn_xonly_any", rec_method = "nn_yonly",
                           options = list(est_cis_window = window,
                                          cand_nn = 20, rec_nn = 1, rec_num_rec = 1,
                                          rec_run_diagnostic = T, est_switch = F,
                                          est_hold_initial = T))

par(mfrow = c(1,2), mar = c(5,5,0.5,0.5))
set.seed(10); plot_umap(dat, xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10); plot_umap(res, multiple_to = "umap_avg", xlab = "UMAP 1", ylab = "UMAP 2")
