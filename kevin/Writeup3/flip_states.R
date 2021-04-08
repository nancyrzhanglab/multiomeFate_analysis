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
dat <- generate_data(obj_next, number_runs = 10, sample_perc = 1, time_tol = 0.01, 
                     verbose = T)
dim(dat$obs_x)

###############################

list_end <- list(which(dat$df_info$time <= 0.1))
vec_start <- which(dat$df_info$time >= 0.9)
set.seed(10)
res <- chromatin_potential(dat$obs_x, dat$obs_y, df_x = dat$df_x, df_y = dat$df_y,
                           vec_start = vec_start, list_end = list_end, 
                           form_method = "average", est_method = "glmnet",
                           cand_method = "nn_xonly_any", rec_method = "nn_yonly",
                           options = list(est_cis_window = window, rec_nn = 10, rec_num_rec = 20))
res$options

par(mfrow = c(1,1))
set.seed(10)
plot_arrow_iteration(res, dat$df_info$time, xlab = "Recruit index",
                     ylab = "True pseudotime")

par(mfrow = c(1,2))
set.seed(10); plot_umap(dat, xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10); plot_umap(res, xlab = "UMAP 1", ylab = "UMAP 2")


## hypothesis: 
# [[we think our algorithm relies a lot of nearest-neighbors,
# so 1) all the arrows would simply be flipped (but still sensible).
# This is in contrast to 2) if all the arrows were now scrambled 
# (some pointing forward, some pointing backward), OR
# 3) the results look exactly the same as beforehand]]
## what we hope would've happened:
# If RNA actually precede ATAC (i.e., the opposite of what 
# we're hypothesis), then we would expect that a hypothetical
# method (that assumes ATAC precedes RNA) to output
# garbage results (i.e., the arrows are pointing in every direction)


# We would want our method to really rely
# on that underlying directionality (ATAC->RNA).

# In SHARE-Seq, the method does not rely on ATAC->RNA but
# their results /support/ this because when they order the cell
# in trajetory (based on ONLY ATAC), 
# they then show that all the arrows (matching ATAC to RNA) 
# point from younger cell to later cells.


# There's actually 3 "aspects" to our problem:
# - In what setting does ATAC procede RNA? : [Specifically, 
## the peaks (that control a gene) open up before said gene expression's
## is high.]
## What we're using is that the ATAC can predict the future RNA.
# - What relation is "constant" for all cells [Currently: we
## are assuming the influence of cis-peaks on the gene is constant]
# - Which cells are the earliest and oldest?

## if we had 1 start state (A) and 2 end states (B and C), and re-run this 
# similar experiment (i.e., we flip the start A with one of the 2 ends, say B)
# a few things could happen:
# - 1) all the arrows go from the "incorrect start" B to the other end state C,
# and none of them go to A
# - 2) all the arrows from the "incorrect start" B to the "incorrect end" A,
# and none of them go to C
# - 3) some of the arrows start from the "incorrect start" B and end up at either A or C.


