rm(list=ls())

set.seed(10)
p1 <- 20; p2 <- 5; genome_length <- 1000; window <- 10
df <- generate_df_simple(p1, p2, genome_length = genome_length, window = window)
mat_g <- generate_gcoef_simple(df$df_x, df$df_y, window = window)
timepoints <- 30

mat_traj <- generate_traj_cascading(df$df_y, timepoints = timepoints)*2
res <- prepare_obj_nextcell(df$df_x, df$df_y, mat_g, list(mat_traj),
                            bool_traj_y = T, verbose = F)

#####
list_traj_mat <- list(mat_traj)
verbose = T
coarseness = 0.1
max_y = 1e5
tol <- 1e-6
mat_x1all <- pmax(pmin(.compute_xfromy(list_traj_mat, mat_g), 1-tol), tol) #x1
mat_y2all <- do.call(rbind, list_traj_mat) #y2

n_total <- nrow(mat_y2all) # count how many unique rows there are
list_time <- list(seq(0, 1, length.out = n_total)) # [note to self: currently hard-coded for linear trajectory]

ht <- hash::hash()
counter <- 1

## grab the relevant rows
## [note to self: hard-code the fact it's a linear trajectory]
i <- 2
idx <- c(max(round(i-coarseness*n_total), 1):min(round(i+coarseness*n_total), n_total-1))
mat_y2 <- mat_y2all[idx,,drop = F] 
mat_x2 <- mat_x1all[idx+1,,drop = F] 

###################

response_prob <- mat_x2
covariate <- mat_y2

x <- covariate
p1 <- ncol(response_prob); p2 <- ncol(covariate)
mat_coef <- matrix(NA, nrow = p2, ncol = p1)
vec_intercept <- rep(NA, length = ncol(response_prob))

i <- 5
if (diff(range(response_prob[,i])) <= tol){
  mat_coef[,i] <- NA; vec_intercept[i] <- median(response_prob[,i])
  
} else {
  y_mat <- cbind(1-response_prob[,i], response_prob[,i])
  fit <- glmnet::glmnet(x = x, y = y_mat, family = "binomial",
                        standardize = FALSE, intercept = TRUE, alpha = 0)
  
  len <- length(fit$lambda)
  mat_coef[,i] <- fit$beta[,len]
  vec_intercept[i] <- fit$a0[len]
}


