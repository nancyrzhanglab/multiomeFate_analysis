rm(list=ls())
set.seed(10)
p1 <- 20; p2 <- 10; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
timepoints <- 20
mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints, max_val = exp(3), min_val = 1)
obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, 
                                 list(mat_traj), verbose = F)
dat <- generate_data(obj_next, number_runs = 10, sample_perc = 0.9, verbose = F)

vec_start <- which(dat$df_info$time <= 0.1)
list_end <- list(which(dat$df_info$time >= 0.9))

####3

set.seed(10)
res <- chromatin_potential(dat$obs_x, dat$obs_y, dat$df_x, dat$df_y,
                           vec_start, list_end, verbose = T)



