rm(list=ls())
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}

p1 <- 100; p2 <- 20; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)

timepoints <- 100
traj <- generate_traj_cascading(df$df_x, timepoints = timepoints)
obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(traj$mat_1), list(traj$mat_2))

################

set.seed(10)
dat <- generate_data(obj_next)

##########################
y <- dat$obs_y[15,]
tmp <- matrix(y, nrow = 1, ncol = length(y))
idx <- RANN::nn2(obj_next$mat_starty, query = tmp, k = 1)$nn.idx[1,1]
idx

image(.rotate(obj_next$mat_starty))
plot(y)
