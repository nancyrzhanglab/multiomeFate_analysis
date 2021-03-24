rm(list=ls())
library(multiomeFate)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

.rotate <- function(mat){t(mat)[,nrow(mat):1]}
 
# step 1: generate the genome
p1 <- 100; p2 <- 20; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
head(df$df_x); head(df$df_y)

#################

# step 2: generate the G matrix, mapping Modality 1 to the next cell's Modality 2
set.seed(10)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window, 
                               signal_sd = 0.1)
dim(mat_g)
mat_g[1:5,1:5]
image(.rotate(mat_g), xlab = "Position in Modality 2", 
      ylab = "Position in Modality 1", asp = T)

#################

# step 3: generate the "blueprint" on how a modality evolves.
# Here, we generate a "cascading from the left" trajectory for Modality 2
# Note: The max_val is exp(max_val) (since these values are the mean, not the natural parameter)
#   and the min_val is 1 (since we will be taking the log of these values in prepare_obj_nextcell)
timepoints <- 100; max_val <- 2
traj_mat <- generate_traj_cascading(df$df_y, timepoints = timepoints, 
                                    max_val = exp(max_val), min_val = 1)
traj_mat <- pmax(traj_mat, 1)
dim(traj_mat)
traj_mat[1:5,1:5]

# we can plot what this trajectory looks like
par(mfrow = c(1,3))
plot(traj_mat[1,], pch = 16, ylim = c(0,exp(max_val)), 
     xlab = "Genomic position", ylab = "Baseline expression", main = "Time: 1")
plot(traj_mat[round(timepoints/2),], pch = 16, ylim = c(0,exp(max_val)), 
     xlab = "Genomic position", ylab = "Baseline expression", main = paste0("Time: ", round(timepoints/2)))
plot(traj_mat[timepoints-1,], pch = 16, ylim = c(0,exp(max_val)), 
     xlab = "Genomic position", ylab = "Baseline expression", main = paste0("Time: ", timepoints))

################

# step 4: compute all the necessary ingredients for the simulation.
# This involves computing the probabilistic (logistic) link from Modality 2 to
# Modality 1.
# This function takes quite a while to compute
set.seed(10)
obj_next <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list_traj_mat = list(traj_mat),
                                 max_y = exp(max_val), verbose = T)

# we can manually test how information from one time point propagates forward.
# This is going under-the-hood. The user never needs to use the following functions.
# to do this, we pick a particular time forward and use the internal-functions to push
#  this point ahead in time. This computes the mean expression of the other modalities/other times.
time <- 40
obj_next$ht[[as.character(time)]]$time
y1_true <- obj_next$mat_y2all[time,]
x1_true <- .bernoulli_xgiveny(y1_true, obj_next$ht[[as.character(time)]]$list_coef$mat_coef, obj_next$ht[[as.character(time)]]$list_coef$vec_intercept)
y2_true <- .possion_ygivenx(x1_true, mat_g)

# plot these values
par(mfrow = c(1,4))
plot(df$df_y$location, y1_true, ylim = c(0, exp(max_val)),  pch = 16,
     xlab = "Genomic position", ylab = "Modal 2 expression", main = "Current time")
plot(df$df_x$location, x1_true, pch = 16,
     xlab = "Genomic position", ylab = "Modal 1 expression", main = "Current time")
plot(df$df_y$location, y2_true, ylim = c(0, exp(max_val)), pch = 16,
     xlab = "Genomic position", ylab = "Modal 1 expression", main = "Next time")
plot(df$df_y$location, y2_true-y1_true, pch = 16,
     xlab = "Genomic position", ylab = "Modal 2 difference", main = "Difference in time")

# we can also generate the random values themselves. 
# This compute the random values in the other modalities/other times.
par(mfrow = c(1,2))
x1 <- .generate_xgiveny(obj_next, y1_true)
plot(x1$x, pch = 16, xlab = "Genomic position", ylab = "Modal 1 expression")
y2 <- .generate_ygivenx(obj_next, x1$x)
plot(y2, pch = 16, xlab = "Genomic position", ylab = "Modal 2 expression")


#################

# step 5: We now generate the data.
set.seed(10)
dat <- generate_data(obj_next, verbose = T)

head(dat$df_info); dim(dat$df_info)
dat$obs_x[1:5,1:5]; dim(dat$obs_x)

# we can make the following plots to see that as pseudotime increases, the total expression increases
par(mfrow = c(1,3))
plot(dat$df_x$location, dat$obs_x[1,], ylim = c(0,1), pch = 16, xlab = "Genomic position",
     ylab = "Modal 1 epxression", main = "Expression at start")
plot(dat$df_info$time, rowSums(dat$obs_x), pch = 16, xlab = "Pseudotime", ylab = "Total expression in Modal 1")
plot(dat$df_info$time, rowSums(dat$obs_y), pch = 16, xlab = "Pseudotime", ylab = "Total expression in Modal 2")

# plot Mode 1
par(mfrow = c(1,2))
plot_activation(dat, mode = 1, reorder = F, cex = 0.2, main = "Ordered by simulation appearance")
plot_activation(dat, mode = 1, reorder = T, cex = 0.2, main = "Ordered by pseudotime")

# plot Mode 2
par(mfrow = c(1,2))
plot_heatmap(dat, mode = 2, reorder = F, main = "Ordered by simulation appearance")
plot_heatmap(dat, mode = 2, reorder = T, main = "Ordered by pseudotime")

# plot UMAP
par(mfrow = c(1,2))
set.seed(10); plot_umap(dat, mode_x = T, mode_y = T, reorder = F, main = "Ordered by simulation appearance")
set.seed(10); plot_umap(dat, mode_x = T, mode_y = T, reorder = T, main = "Ordered by pseudotime")
