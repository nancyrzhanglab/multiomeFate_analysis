rm(list=ls())
load("tests/assets/test.RData")

bandwidth <- 200
discretization_stepsize <- 10

frag_win <- .extract_fragment_from_cutmat(cutmat_winning)
frag_die <- .extract_fragment_from_cutmat(cutmat_dying)
frag_all <- c(frag_win, frag_die)

trials <- 10

i <- 1
set.seed(10*i)
idx <- sample(1:length(frag_all), size = round(length(frag_all)/2))
idx_die <- sample(1:length(idx), size = round(length(idx)/2))
idx_win <- sample(1:(length(frag_all) - length(idx)),
                  size = round((length(frag_all) - length(idx))/2))

# res <- .compute_crossfit_teststat(
#   bandwidth = bandwidth,
#   frag_die = frag_all[idx],
#   frag_win = frag_all[-idx],
#   idx_die = idx_die,
#   idx_win = idx_win,
#   peak_locations = peak_locations,
#   peak_prior = peak_prior,
#   peak_width = peak_width,
#   discretization_stepsize = discretization_stepsize,
#   bool_lock_within_peak = T,
#   max_iter = 100,
#   min_fragments = 6,
#   min_prior = 0.01,
#   num_peak_limit = 4,
#   tol = 1e-6,
#   verbose = 0
# )
# res$teststat

############################

bool_lock_within_peak = T 
max_iter = 100
min_fragments = 6
min_prior = 0.01
num_peak_limit = 4
tol = 1e-6
verbose = 0
frag_die = frag_all[idx] 
frag_win = frag_all[-idx]

# fit1 <- .lrt_onefold(
#   bandwidth = bandwidth,
#   frag_die = frag_die, 
#   frag_win = frag_win,
#   idx_die_split1 = idx_die,
#   idx_win_split1 = idx_win,
#   peak_locations = peak_locations,
#   peak_prior = peak_prior,
#   peak_width = peak_width,
#   discretization_stepsize = discretization_stepsize, 
#   bool_lock_within_peak = bool_lock_within_peak, 
#   max_iter = max_iter,
#   min_prior = min_prior,
#   num_peak_limit = num_peak_limit,
#   tol = tol,
#   verbose = verbose
# )

#################
idx_die_split1 = idx_die
idx_win_split1 = idx_win
len_die <- length(frag_die); len_win <- length(frag_win)

fit_win <- peak_mixture_modeling(
  bandwidth = bandwidth,
  cutmat = NULL, 
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  peak_width = peak_width,
  discretization_stepsize = discretization_stepsize, 
  bool_lock_within_peak = bool_lock_within_peak, 
  bool_freeze_prior = T, # assumed to not have a lot of fragments, so we need to freeze
  fragment_locations = frag_win[idx_win_split1], 
  max_iter = max_iter,
  min_prior = min_prior,
  num_peak_limit = num_peak_limit,
  return_dist_mat = F,
  tol = tol,
  verbose = verbose
)

print("fit die")
fit_die <- peak_mixture_modeling(
  bandwidth = bandwidth,
  cutmat = NULL, 
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  peak_width = peak_width,
  discretization_stepsize = discretization_stepsize, 
  bool_lock_within_peak = bool_lock_within_peak, 
  bool_freeze_prior = T,
  fragment_locations = frag_die[idx_die_split1], 
  max_iter = max_iter,
  min_prior = min_prior,
  num_peak_limit = num_peak_limit,
  return_dist_mat = F,
  tol = tol,
  verbose = verbose
)

# p0 is for the null
# compute the win=die on the second fold
print("fit both")
fit_both <- peak_mixture_modeling(
  bandwidth = bandwidth,
  cutmat = NULL, 
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  peak_width = peak_width,
  discretization_stepsize = discretization_stepsize, 
  bool_lock_within_peak = bool_lock_within_peak, 
  bool_freeze_prior = T,
  fragment_locations = c(frag_win[-idx_win_split1], frag_die[-idx_die_split1]), 
  max_iter = max_iter,
  min_prior = min_prior,
  num_peak_limit = num_peak_limit,
  return_dist_mat = T, # needed for the numerator likelihood later
  tol = tol,
  verbose = verbose
)
stopifnot(nrow(fit_both$dist_mat) == len_win + len_die - length(idx_win_split1) - length(idx_die_split1))
loglikelihood_denom <- fit_both$loglikelihood_val

# compute the likelihood ratio of L(win!=die)/L(win=die) on the second fold
# that is, p1(Y_second)/p0(Y_second)
loglikelihood_outofsample_win <- .compute_loglikelihood(
  dist_mat = fit_both$dist_mat[1:(len_win - length(idx_win_split1)),],
  grenander_obj = fit_win$grenander_obj,
  prior_vec = fit_win$prior_vec
)
loglikelihood_outofsample_die <- .compute_loglikelihood(
  dist_mat = fit_both$dist_mat[(len_win-length(idx_win_split1)+1):nrow(fit_both$dist_mat),],
  grenander_obj = fit_die$grenander_obj,
  prior_vec = fit_die$prior_vec
)
loglikelihood_num <- loglikelihood_outofsample_win + loglikelihood_outofsample_die

# return test statistic
teststat <- exp(loglikelihood_num - loglikelihood_denom)

#####################

plot(fit_win$grenander_obj$x, fit_win$grenander_obj$pdf, main = "Win")
plot(fit_die$grenander_obj$x, fit_die$grenander_obj$pdf, main = "Die")
plot(fit_both$grenander_obj$x, fit_both$grenander_obj$pdf, main = "Both")
