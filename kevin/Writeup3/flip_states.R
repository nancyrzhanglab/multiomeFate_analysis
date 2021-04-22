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


###############################

list_end <- list(which(dat$df_info$time <= 0.1))
vec_start <- which(dat$df_info$time >= 0.9)
set.seed(10)
# run the estimator
res <- chromatin_potential(dat$obs_x, dat$obs_y, df_x = dat$df_x, df_y = dat$df_y,
                           vec_start = vec_start, list_end = list_end, 
                           form_method = "average", est_method = "glmnet",
                           cand_method = "nn_xonly_avg", rec_method = "nn_yonly",
                           options = list(est_cis_window = window, cand_num_cand = 30,
                                          cand_nn = 30,
                                          rec_nn = 10, rec_num_rec = 10,
                                          rec_run_diagnostic = T))
res$options

par(mfrow = c(1,1), mar = c(5,5,0.5,0.5))
set.seed(10)
plot_arrow_iteration(res, dat$df_info$time, xlab = "Recruit index",
                     ylab = "True pseudotime")

par(mfrow = c(1,2), mar = c(5,5,0.5,0.5))
set.seed(10); plot_umap(dat, ghost_neighbor = res$ht_neighbor, xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10); plot_umap(res, multiple_to = "ghost", xlab = "UMAP 1", ylab = "UMAP 2")

par(mfrow = c(1,2), mar = c(5,5,3,0.5))
set.seed(10); plot_umap(res, multiple_to = "ghost", xlab = "UMAP 1", ylab = "UMAP 2", 
                        percent_arrows = 1, num_col_arrows = 10, col_arrows_by = "order_rec",
                        main = "Order of recruitment")
set.seed(10); plot_umap(res, multiple_to = "ghost", xlab = "UMAP 1", ylab = "UMAP 2", 
                        percent_arrows = 1, col_arrows_by = "direction", vec_time = dat$df_info$time,
                        main = "Pseudotime direction")

par(mfrow = c(1,1))
plot(jitter(res$df_res$order_rec), jitter(res$df_res$num_cand), pch = 16, 
     xlab = "Recruitment index", ylab = "# times considered as candidate",
     col = rgb(0.5,0.5,0.5,0.5))


par(mfrow = c(1,2), mar = c(5,5,0.5,0.5))
set.seed(10); plot_umap(dat, xlab = "UMAP 1", ylab = "UMAP 2")
set.seed(10); plot_umap(res, multiple_to = "umap_avg", xlab = "UMAP 1", ylab = "UMAP 2")

par(mfrow = c(1,2), mar = c(5,5,3,0.5))
set.seed(10); plot_umap(res, multiple_to = "umap_avg", xlab = "UMAP 1", ylab = "UMAP 2", 
                        percent_arrows = 1, num_col_arrows = 10, col_arrows_by = "order_rec",
                        main = "Order of recruitment")
set.seed(10); plot_umap(res, multiple_to = "umap_avg", xlab = "UMAP 1", ylab = "UMAP 2", 
                        percent_arrows = 1, col_arrows_by = "direction", vec_time = dat$df_info$time,
                        main = "Pseudotime direction")

####################
length(res$list_diagnos)
for(iter in rev(1:length(res$list_diagnos))){
  print(iter)
  tmp <- res$list_diagnos[[as.character(iter)]]$recruit$postprocess
  tmp <- t(tmp[,2:4])
  tmp <- as.table(tmp)
  colnames(tmp) <- res$list_diagnos[[iter]]$recruit$postprocess[,1]
  barplot(tmp, col=c("darkblue", "red", "green"), horiz = F,
          main = paste0("Flipped: Iteration ", iter))
  idx <- which(res$list_diagnos[[as.character(iter)]]$recruit$postprocess$selected)
  for(i in idx){
    start <- 0.2+1.2*(i-1)
    mid <- start+0.5
    lines(rep(mid,2), c(-1e4,1e4), col = 'purple', lty = 2, lwd = 2)
  }
}

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


