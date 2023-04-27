rm(list=ls())
set.seed(10)
peak_locations <- c(1500, 1900)
names(peak_locations) <- paste0("p:", 1:2)
num_frags <- 120
peak_width <- 100
peak_prior <- c(0.5,0.5)

frag_location <- rep(c(1400,1450,1850), each = num_frags/3)
cutmat_dying <- sapply(frag_location, function(x){
  vec <- rep(0, length = 1000)
  names(vec) <- c(1001:2000)
  vec[as.character(x)] <- 1
  vec
})
cutmat_dying <- Matrix::Matrix(t(cutmat_dying), sparse = T)
rownames(cutmat_dying) <- paste0("cell:", 1:nrow(cutmat_dying))

frag_location <- rep(c(1100,1250,1600), each = num_frags/3)
cutmat_winning <- sapply(frag_location, function(x){
  vec <- rep(0, length = 1000)
  names(vec) <- c(1001:2000)
  vec[as.character(x)] <- 1
  vec
})
cutmat_winning <- Matrix::Matrix(t(cutmat_winning), sparse = T)
rownames(cutmat_winning) <- paste0("cell:", 1:nrow(cutmat_winning))

res <- peak_testing(
  cutmat_dying = cutmat_dying, 
  cutmat_winning = cutmat_winning,
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  peak_width = peak_width,
  verbose = 1
)

res$teststat
res$pvalue
