rm(list=ls())
load("tests/assets/test.RData")
res <- peak_mixture_modeling(bin_limits = bin_limits,
                             bin_midpoints = bin_midpoints, 
                             cutmat = cutmat_dying, 
                             peak_locations = peak_locations,
                             peak_prior = peak_prior,
                             bool_freeze_prior = F,
                             return_assignment_intial = T,
                             return_assignment_mat = T,
                             return_bin_mat = T,
                             return_theta_initial = T,
                             verbose = 3)
expect_true(inherits(res, "peakDistribution"))

round(res$theta_vec, 2)
head(res$bin_mat)
head(round(res$assignment_mat,2))
round(res$prior_vec, 2)

idx <- which.max(res$prior_vec)
res$bin_mat[,idx]
plot(res$bin_mat[,idx])
for(i in 1:length(res$prior_vec)){
  par(mfrow = c(1,2))
  plot(jitter(res$bin_mat[,i]), res$assignment_mat[,i], 
       pch = 16, main = paste0("Fit: ",i, ", Prior: ", round(res$prior_vec[i],2)))
  plot(jitter(res$bin_mat[,i]), res$assignment_init[,i], 
       pch = 16, main = paste0("Initial: ",i, ", Prior: ", round(peak_prior[i],2)))
}
