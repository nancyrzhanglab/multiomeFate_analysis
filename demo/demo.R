rm(list=ls())
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}
 
p1 <- 100; p2 <- 20; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
head(df$df_x); head(df$df_y)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
dim(mat_g)
mat_g[1:5,1:5]
image(.rotate(mat_g))

timepoints <- 100
traj <- generate_traj_cascading(df$df_x, timepoints = timepoints)
names(traj); dim(traj$mat_1)
traj$mat_1[1:5,1:5]
traj$mat_2[1:5,1:5]
par(mfrow = c(1,3))
plot(traj$mat_1[1,], pch = 16, ylim = c(0,1))
plot(traj$mat_1[round(timepoints/2),], pch = 16, ylim = c(0,1))
plot(traj$mat_1[timepoints-1,], pch = 16, ylim = c(0,1))

################
obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(traj$mat_1), list(traj$mat_2))

time <- 40
obj_next$ht[[as.character(time)]]$time
x1_true <- traj$mat_1[time,]
y2_true <- .possion_ygivenx(x1_true, mat_g)
x2_true <- .bernoulli_xgiveny(y2_true, obj_next$ht[[as.character(time)]]$list_coef$mat_coef, obj_next$ht[[as.character(time)]]$list_coef$vec_intercept)

par(mfrow = c(1,4))
plot(df$df_x$location, x1_true, ylim = c(0,1))
plot(df$df_y$location, y2_true)
plot(df$df_x$location, x2_true, ylim = c(0,1))
plot(df$df_x$location, x2_true-x1_true)

par(mfrow = c(1,2))
y2 <- generate_ygivenx(obj_next, x1_true)
plot(y2)
x2 <- generate_xgiveny(obj_next, y2)
plot(x2$x)

#################

set.seed(10)
dat <- generate_data(obj_next)

head(dat$df_info); dim(dat$df_info)
dat$obs_x[1:5,1:5]; dim(dat$obs_x)
plot(dat$df_x$location, dat$obs_x[1,], ylim = c(0,1))
plot(dat$df_info$counter, rowSums(dat$obs_x))
plot(dat$df_info$counter, rowSums(dat$obs_y))

par(mfrow = c(1,2))
plot_activation(dat, reorder = F, cex = 0.2)
plot_activation(dat, reorder = T, cex = 0.2)

par(mfrow = c(1,2))
plot_heatmap(dat, reorder = F)
plot_heatmap(dat, reorder = T)


par(mfrow = c(1,2))
set.seed(10); plot_umap(dat, color_by = "counter")
set.seed(10); plot_umap(dat, color_by = "time")
