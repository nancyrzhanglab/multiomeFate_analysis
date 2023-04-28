rm(list=ls())
load("tests/assets/test.RData")
idx <- which(peak_prior <= 0.05)
peak_locations <- peak_locations[-idx]
peak_prior <- peak_prior[-idx]
peak_prior <- peak_prior/sum(peak_prior)

Rprof("profile.out")
time_start <- Sys.time()
# Call the function to be profiled
res_win <- peak_mixture_modeling(cutmat = cutmat_dying,
                                 peak_locations = peak_locations,
                                 peak_prior = peak_prior,
                                 peak_width = peak_width,
                                 return_lowerbound = T,
                                 return_dist_mat = T,
                                 max_iter = 100,
                                 verbose = 0)
res_die <- peak_mixture_modeling(cutmat = cutmat_winning,
                                 peak_locations = peak_locations,
                                 peak_prior = peak_prior,
                                 peak_width = peak_width,
                                 return_dist_mat = T,
                                 max_iter = 100,
                                 verbose = 0)
res_both <- peak_mixture_modeling(cutmat = rbind(cutmat_winning,
                                                 cutmat_dying),
                                  peak_locations = peak_locations,
                                  peak_prior = peak_prior,
                                  peak_width = peak_width,
                                  return_dist_mat = T,
                                  max_iter = 100,
                                  verbose = 0)
time_end <- Sys.time()
Rprof(NULL)
time_end - time_start

zz <- summaryRprof("profile.out")
head(zz$by.self, 20)
head(zz$by.total, 20)

#################

Rprof("profile2.out")
time_start2 <- Sys.time()
# Call the function to be profiled
res_win <- peak_mixture_modeling(cutmat = cutmat_dying,
                                 peak_locations = peak_locations,
                                 peak_prior = peak_prior,
                                 peak_width = peak_width,
                                 return_dist_mat = T,
                                 return_lowerbound = F,
                                 max_iter = 100,
                                 verbose = 0)
res_die <- peak_mixture_modeling(cutmat = cutmat_winning,
                                 peak_locations = peak_locations,
                                 peak_prior = peak_prior,
                                 peak_width = peak_width,
                                 return_dist_mat = T,
                                 return_lowerbound = F,
                                 max_iter = 100,
                                 verbose = 0)
res_both <- peak_mixture_modeling(cutmat = rbind(cutmat_winning,
                                                 cutmat_dying),
                                  peak_locations = peak_locations,
                                  peak_prior = peak_prior,
                                  peak_width = peak_width,
                                  return_dist_mat = T,
                                  return_lowerbound = F,
                                  max_iter = 100,
                                  verbose = 0)
time_end2 <- Sys.time()
Rprof(NULL)
time_end2 - time_start2
