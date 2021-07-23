rm(list=ls()); gc(T)
load("../../../../out/kevin/Writeup3b/10x_embryo.RData")

vec_start <- which(celltype == "Radial glia")
list_end <- list(which(celltype == "Oligodendrocyte"),
                 which(celltype == "Forebrain GABAergic"),
                 which(celltype == "Cortical or hippocampal glutamatergic"))
(length(vec_start) + length(unlist(list_end)))/nrow(mat_x)

rank_x <- 30
rank_y <- 50
df_x <- data.frame(name = colnames(mat_x))
df_y <- data.frame(name = colnames(mat_y))
mat_x <- as.matrix(mat_x)
mat_y <- as.matrix(mat_y)

set.seed(10)
prep_obj <- multiomeFate::chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                                      vec_start, list_end,
                                                      form_method = "average_weighted",
                                                      est_method = "threshold_glmnet",
                                                      rec_method = "distant_cor",
                                                      ht_map = ht_map,
                                                      options = list(nn_nn = 30, dim_nlatent_x = rank_x,
                                                                     dim_nlatent_y = rank_y, 
                                                                     est_num_iterations = 4,
                                                                     rec_bool_pred_nn = T,
                                                                     est_cv_choice = "lambda.min",
                                                                     form_stepsize = 0.5,
                                                                     form_min_weight = 0,
                                                                     est_verbose = T))


#########

options <- prep_obj$options
form_options <- options$form_options
est_options <- options$est_options
df_res <- prep_obj$df_res
tmp <- multiomeFate:::.init_est_matrices(mat_x, mat_y, df_res, form_options)
mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
weights <-  multiomeFate:::.init_weights(nrow(mat_x1), form_options)
p1 <- ncol(mat_x1); p2 <- ncol(mat_y2)

j <- 527
idx_x <- est_options$ht_map[[as.character(j)]]
idx_y <- j
# multiomeFate:::.threshold_glmnet_fancy(mat_x1[,idx_x,drop = F], mat_y2[,idx_y],
#                         weights = weights,
#                         family = est_options$family, 
#                         switch = est_options$switch, switch_cutoff = est_options$switch_cutoff,
#                         alpha = est_options$alpha, intercept = est_options$intercept,
#                         cv = est_options$cv, nfolds = est_options$nfolds, cv_choice = est_options$cv_choice,
#                         num_iterations = est_options$num_iterations, 
#                         initial_quantile = est_options$initial_quantile)

#################
x <- mat_x1[,idx_x,drop = F]
y <- mat_y2[,idx_y]
family = est_options$family
switch = est_options$switch
switch_cutoff = est_options$switch_cutoff
alpha = est_options$alpha
intercept = est_options$intercept
cv = est_options$cv
nfolds = est_options$nfolds
cv_choice = est_options$cv_choice
num_iterations = est_options$num_iterations
initial_quantile = est_options$initial_quantile

prev_threshold <- stats::quantile(y, probs = initial_quantile)
iter <- 1

while(iter <= num_iterations){
  # update the regression
  print("asdf1")
  idx <- which(y > prev_threshold)
  if(length(idx) > 0){
    if(length(weights) == 1) weight_vec <- NA else weight_vec <- weights[idx]
    res_glm <- multiomeFate:::.glmnet_fancy(x[idx,,drop = F], y[idx], weight_vec,
                             family, switch, switch_cutoff,
                             alpha, intercept,
                             cv, nfolds, cv_choice)
  } else {
    return(list(val_int = 0, vec_coef = rep(0, ncol(x)),
                val_threshold = prev_threshold))
  }
  
  
  # update the threshold
  print("asdf2")
  next_threshold <- multiomeFate:::.update_threshold_glmnet(x, y, weights, res_glm)
  
  if(abs(next_threshold - prev_threshold) <= tol) break()
  iter <- iter+1
}

##################

x <- mat_x; y <- vec_y
n <- length(y); p <- ncol(x)
weights <- weight_vec
res <- glmnet::cv.glmnet(x, y, weights = weights, family = family, nfolds = nfolds, alpha = alpha,
                         intercept = intercept)
res <- glmnet::cv.glmnet(x, y)
res <- glmnet::glmnet(x, y)



