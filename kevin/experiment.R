rm(list=ls())
load("tests/assets/test.RData")
bandwidth <- 200
discretization_stepsize <- 10
cutmat_all <- rbind(cutmat_dying, cutmat_winning)
trials <- 10

i <- 1
set.seed(10*i)
idx1 <- sample(1:nrow(cutmat_all), round(nrow(cutmat_all)/2))

res <- peak_testing(
  bandwidth = bandwidth,
  cutmat_dying = cutmat_all[idx1,],
  cutmat_winning = cutmat_all[-idx1,],
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  peak_width = peak_width,
  discretization_stepsize = discretization_stepsize,
  verbose = 0
)

###########################

cutmat_dying = cutmat_all[idx1,]
cutmat_winning = cutmat_all[-idx1,]
bool_lock_within_peak = T
max_iter = 100
min_fragments = 6
min_prior = 0.01
num_peak_limit = 4
tol = 1e-6
verbose = 1

frag_win <- .extract_fragment_from_cutmat(cutmat_winning)
frag_die <- .extract_fragment_from_cutmat(cutmat_dying)

idx_win <- sample(1:length(frag_win), size = round(length(frag_win)/2))
idx_die <- sample(1:length(frag_die), size = round(length(frag_die)/2))

idx_die_split1 = idx_die
idx_win_split1 = idx_win

# fit_die <- peak_mixture_modeling(
#   bandwidth = bandwidth,
#   cutmat = NULL, 
#   peak_locations = peak_locations,
#   peak_prior = peak_prior,
#   peak_width = peak_width,
#   discretization_stepsize = discretization_stepsize, 
#   bool_lock_within_peak = bool_lock_within_peak, 
#   bool_freeze_prior = T,
#   fragment_locations = frag_die[idx_die_split1], 
#   max_iter = max_iter,
#   min_prior = min_prior,
#   num_peak_limit = num_peak_limit,
#   return_dist_mat = F,
#   tol = tol,
#   verbose = verbose
# )

##################

fragment_locations = frag_die[idx_die_split1]
cutmat <- NULL
return_assignment_mat = F
return_dist_mat = F
bool_freeze_prior = T

dist_mat <- .compute_frag_peak_matrix(
  bool_lock_within_peak = bool_lock_within_peak,
  cutmat = cutmat,
  fragment_locations = fragment_locations,
  num_peak_limit = num_peak_limit,
  peak_locations = peak_locations,
  peak_width = peak_width
)
num_frags <- nrow(dist_mat)
if(is.na(discretization_stepsize)) discretization_stepsize <- max(round(max(dist_mat@x)/2000),5)

grenander_obj <- .initialize_grenander(bandwidth = bandwidth,
                                       dist_mat = dist_mat,
                                       discretization_stepsize = discretization_stepsize)

iter <- 1
loglikelihood_vec <- .compute_loglikelihood(
  dist_mat = dist_mat,
  grenander_obj = grenander_obj,
  prior_vec = peak_prior
)
if(verbose > 1) print(paste0("Initial log-likelihood: ", round(loglikelihood_vec[1],2)))
prior_vec <- peak_prior

# TODO: Return if there are no fragments

while(TRUE){
  if(iter > max_iter) break()
  
  if(verbose > 1) print("E-step")
  assignment_mat <- .e_step(
    dist_mat = dist_mat,
    grenander_obj = grenander_obj,
    prior_vec = prior_vec
  )
  print(length(assignment_mat@x))
  
  if(verbose > 1) print("M-step")
  grenander_obj_new <- .m_step(
    assignment_mat = assignment_mat,
    bandwidth = bandwidth,
    dist_mat = dist_mat,
    discretization_stepsize = discretization_stepsize
  )
  if(!bool_freeze_prior) { prior_vec <- .compute_prior(assignment_mat = assignment_mat, min_prior = min_prior) }
  
  if(verbose) print("Computing likelihood")
  loglikelihood_val <- .compute_loglikelihood(
    dist_mat = dist_mat,
    grenander_obj = grenander_obj_new,
    prior_vec = peak_prior
  )
  
  loglikelihood_vec <- c(loglikelihood_vec, loglikelihood_val)
  iter <- length(loglikelihood_vec)
  if(length(loglikelihood_vec) >= 2){
    if(verbose > 0) print(paste0("Iteration: ", iter, ", log-likelihood: ", round(loglikelihood_vec[iter],2)))
    if(abs(loglikelihood_vec[iter] - loglikelihood_vec[iter-1]) <= tol) break()
  }
  grenander_obj <- grenander_obj_new
}

###########

# let's see where it crashed
grenander_obj_new <- .m_step(
  assignment_mat = assignment_mat,
  bandwidth = bandwidth,
  dist_mat = dist_mat,
  discretization_stepsize = discretization_stepsize
)
length(assignment_mat@x)
length(dist_mat@x)
