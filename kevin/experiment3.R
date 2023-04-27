rm(list=ls())
load("tests/assets/test.RData")
idx <- which(peak_prior <= 0.05)
peak_locations <- peak_locations[-idx]
peak_prior <- peak_prior[-idx]
peak_prior <- peak_prior/sum(peak_prior)

res_win <- peak_mixture_modeling(cutmat = cutmat_dying,
                                 peak_locations = peak_locations,
                                 peak_prior = peak_prior,
                                 peak_width = peak_width,
                                 return_dist_mat = T,
                                 max_iter = 100,
                                 verbose = 1)
res_die <- peak_mixture_modeling(cutmat = cutmat_winning,
                                 peak_locations = peak_locations,
                                 peak_prior = peak_prior,
                                 peak_width = peak_width,
                                 return_dist_mat = T,
                                 max_iter = 100,
                                 verbose = 1)

res_both <- peak_mixture_modeling(cutmat = rbind(cutmat_winning,
                                                 cutmat_dying),
                                  peak_locations = peak_locations,
                                  peak_prior = peak_prior,
                                  peak_width = peak_width,
                                  return_dist_mat = T,
                                  max_iter = 100,
                                  verbose = 1)

lrt <- -2 * (res_both$loglikelihood_val - (res_win$loglikelihood_val + res_die$loglikelihood_val))
pval <- 1-stats::pchisq(lrt, df = 100)

#####################

x_max <- min(max(res_win$grenander_obj$x), max(res_die$grenander_obj$x))
bin_cutoff <- exp(seq(log(0.1),log(x_max-1), length.out = 6))[-1]
res_win_coarsen <- coarsen_density(res_win$grenander_obj, bin_cutoff = bin_cutoff)
res_die_coarsen <- coarsen_density(res_die$grenander_obj, bin_cutoff = bin_cutoff)
res_both_coarsen <- coarsen_density(res_both$grenander_obj, bin_cutoff = bin_cutoff)

plot(res_win_coarsen$x, res_win_coarsen$pdf)
plot(res_die_coarsen$x, res_die_coarsen$pdf)
plot(res_both_coarsen$x, res_both_coarsen$pdf)

ll_win <- .compute_loglikelihood(dist_mat = res_win$dist_mat,
                                 grenander_obj =  res_win_coarsen, 
                                 log_prior_vec = log(res_win$prior_vec))
ll_die <- .compute_loglikelihood(dist_mat = res_die$dist_mat,
                                 grenander_obj =  res_die_coarsen, 
                                 log_prior_vec = log(res_die$prior_vec))
ll_both <- .compute_loglikelihood(dist_mat = res_both$dist_mat,
                                  grenander_obj =  res_both_coarsen, 
                                  log_prior_vec = log(res_both$prior_vec))

lrt <- -2 * (ll_both - (ll_win + ll_die))
pval <- 1-stats::pchisq(lrt, df = length(res_win_coarsen$pdf)-1)

list(lrt = lrt, pval = pval)
