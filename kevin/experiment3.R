set.seed(10)
peak_locations <- c(1500, 1900)
names(peak_locations) <- paste0("p:", 1:2)
num_frags <- 120
frag_location <- rep(c(1100,1300,1800), each = num_frags/3)
cutmat <- sapply(frag_location, function(x){
  vec <- rep(0, length = 1000)
  names(vec) <- c(1001:2000)
  vec[as.character(x)] <- 1
  vec
})
cutmat <- Matrix::Matrix(t(cutmat), sparse = T)
rownames(cutmat) <- paste0("cell:", 1:nrow(cutmat))
peak_prior <- c(0.5, 0.5)
peak_width <- 250

res <- peak_mixture_modeling(cutmat = cutmat,
                             peak_locations = peak_locations,
                             peak_prior = peak_prior,
                             peak_width = peak_width,
                             bool_freeze_prior = F,
                             min_prior = 0,
                             return_assignment_mat = T,
                             return_dist_mat = T,
                             return_lowerbound = T,
                             verbose = 0)

res$assignment_mat
