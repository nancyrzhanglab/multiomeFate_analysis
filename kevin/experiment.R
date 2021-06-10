rm(list=ls())
set.seed(10)
g <- igraph::graph_from_edgelist(matrix(c(4,1, 4,5, 2,5, 3,5), nrow = 4, ncol = 2, byrow = T), 
                                 directed = F)
g <- igraph::set_vertex_attr(g, name = "lag", index = 4, value = 3)
g <- igraph::set_vertex_attr(g, name = "lag", index = 5, value = 5)
idx_root <- 4; num_waves <- 10; num_per_wave <- 5; distinct_waves <- 2
combn_wave_mat <- simulate_combn_wave_mat(g, idx_root, num_waves = num_waves,
                                          num_per_wave = num_per_wave, 
                                          distinct_waves = distinct_waves)

res <- simulate_data_input(combn_wave_mat, 
                           x_exp_baseline = 0.1, x_exp_max = 0.7,
                           x_sd_biological = 0.1, x_sd_technical = 0.5, 
                           y_exp_baseline = 0.1, y_sd_technical = 2,
                           num_unrelated_x = 100, num_unrelated_y = 50, 
                           time_on = 15, time_windup = 15, 
                           max_lag = 20, min_lag = 10,
                           x_unrelated_intervals = 2,
                           x_unrelated_max = 0.1, y_unrelated_max = 2)
df_x <- res$df_x; df_y <- res$df_y
list_xnoise <- res$list_xnoise; list_ynoise <- res$list_ynoise

df_cell <- simulate_df_cell(1000, time_max = max(df_y$time_end_scaffold, na.rm = T),
                            num_branch = 3)

set.seed(10)
dat <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)
idx <- which(df_cell$time >= 90)
dat$obs_x <- dat$obs_x[-idx,]; dat$obs_y <- dat$obs_y[-idx,]
dat$mean_x <- dat$mean_x[-idx,]; dat$mean_y <- dat$mean_y[-idx,]
df_cell <- df_cell[-idx,]

mat_x <- dat$obs_x; mat_y <- dat$obs_y

##########################3

vec_start <- which(df_cell$time <= 10)
list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
  intersect(which(df_cell$branch == branch), which(df_cell$time >= 80))
})
# zz <- svd(mat_x); plot(zz$d[2:50])
rank_x <- 15
# zz <- svd(mat_y); plot(zz$d[2:50])
rank_y <- 15
set.seed(10)
prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                        vec_start, list_end,
                                        form_method = "average",
                                        est_method = "threshold_glmnet",
                                        rec_method = "distant_cor",
                                        options = list(nn_nn = 20, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, est_cis_window = 30,
                                                       est_num_iterations = 4,
                                                       rec_bool_pred_nn = T))

set.seed(10)
# res <- chromatin_potential(prep_obj,  verbose = T, bool_oracle = F)

#####################
df_cell2 = df_cell
df_cell = NA
verbose = T
bool_oracle = F
mat_g_init = NA
vec_g_init = rep(0, ncol(mat_y))
vec_threshold_init = rep(0, ncol(mat_y))

# pull the appropriate objects for convenience
mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
df_x <- prep_obj$df_x; df_y <- prep_obj$df_y
df_res <- prep_obj$df_res; dim_reduc_obj <- prep_obj$dim_reduc_obj
nn_g <- prep_obj$nn_g; nn_mat <- prep_obj$nn_mat; nn_obj <- prep_obj$nn_obj
list_diagnos <- prep_obj$list_diagnos; options <- prep_obj$options

dim_options <- options$dim_options; nn_options <- options$nn_options
form_options <- options$form_options; est_options <- options$est_options
cand_options <- options$cand_options; rec_options <- options$rec_options

# initialize
n <- nrow(mat_x)

ht_neighbor <- .init_chrom_ht(which(df_res$order_rec == 0))
tmp <- .init_est_matrices(mat_x, mat_y, df_res, form_options)
mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
list_diagnos <- list()
iter <- 1

if(verbose) print(paste0("Iteration ", iter, ": Recruited percentage (", 
                         round(sum(!is.na(df_res$order_rec))/nrow(df_res), 2), ")"))
## estimate res_g
if((iter == 1 | est_options$hold_initial) && !any(is.na(mat_g_init)) && !any(is.na(vec_g_init))){
  res_g <- list(mat_g = mat_g_init, vec_g = vec_g_init, vec_threshold = vec_threshold_init)
} else {
  res_g <- .estimate_g(mat_x1, mat_y2, est_options)
}

## construct candidate set
res_cand <- .candidate_set(mat_x, mat_y, df_res, nn_mat, cand_options)
df_res <- .update_chrom_df_cand(df_res, res_cand$vec_cand)
stopifnot(all(is.na(df_res$order_rec[res_cand$vec_cand])))
list_diagnos[[as.character(iter)]]$candidate <- res_cand$diagnostic

## recruit an element from the candidate se
enforce_matched <- length(which(df_res$order_rec == 0)) > length(which(df_res$order_rec > 0)) & !bool_oracle
# res_rec <- .recruit_next(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
#                          dim_reduc_obj, nn_g, nn_mat, nn_obj, enforce_matched,
#                          df_cell, rec_options)
# 
# for(i in 1:length(res_rec$rec)){
#   print(table(df_cell2$branch[res_rec$rec[[i]]$from]))
#   print(table(df_cell2$branch[res_rec$rec[[i]]$to]))
#   print("===")
# }

##############3
vec_cand = res_cand$vec_cand
nn_size <- ncol(nn_mat)

if(enforce_matched){
  matched_idx <- which(!is.na(df_res$order_rec))
  table(df_cell2$branch[matched_idx])
} else {
  matched_idx <- NA
}

# apply mat_g to mat_x
if(rec_options$bool_avg_from){
  list_nn <- .extract_nn_list(vec_cand, nn_mat)
  mat_avg_x <- .construct_avg_mat(mat_x, list_nn)
  mat_avg_y <- .construct_avg_mat(mat_y, list_nn)
  pred_y <- .predict_yfromx(mat_avg_x, res_g, rec_options$family)
} else {
  list_nn <- lapply(vec_cand, function(x){x})
  mat_avg_x <- mat_x[vec_cand,]
  mat_avg_y <- mat_y[vec_cand,]
  pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g, rec_options$family)
}

# initialize variables for the loop
if(!rec_options$parallel && future::nbrOfWorkers() == 1){
  my_lapply <- pbapply::pblapply
  if(rec_options$verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
} else {
  my_lapply <- future.apply::future_lapply
}

i=1
cell <- vec_cand[i]
vec <- c(.apply_dimred(mat_avg_x[i,], dim_reduc_obj$x),
         .apply_dimred(pred_y[i,], dim_reduc_obj$y))

# allow cell to be matched to potentially any other cell
if(!rec_options$bool_pred_nn){
  if(!enforce_matched){
    nn_pred <- nn_obj$getNNsByVector(vec, nn_size) + 1
  } else {
    nn_pred <- sample(matched_idx, size = ceiling(rec_options$matched_sampling_rate * length(matched_idx)))
  }
  list_nn_to <- lapply(nn_pred, function(x){x})
  
  # match to only cells near target cell
} else {
  list_nn_to <- .find_to_list_matched(mat_x, mat_y, cell, dim_reduc_obj, 
                                      include_idx = setdiff(matched_idx, list_nn[[i]]), 
                                      exclude_idx = list_nn[[i]], nn_mat,
                                      rec_options)
}

zz <- sapply(list_nn_to, function(x){
  tmp <- table(df_cell2$branch[x])
  as.numeric(names(tmp)[1])
})
df_cell2$branch[cell]

plot(mat_y[cell,], col = as.numeric(as.factor(df_y$branch)), pch = 16)
plot(pred_y[i,], col = as.numeric(as.factor(df_y$branch)), pch = 16)
plot(mat_x[i,], col = as.numeric(as.factor(df_x$branch)), pch = 16)
plot(mat_y[266,], col = as.numeric(as.factor(df_y$branch)), pch = 16)

# from this set of cells, find the ones with highest pearson
pred_diff <- pred_y[i,] - mat_avg_y[i,]
mat_nn_to <- .construct_avg_mat(mat_y, list_nn_to)
cor_vec <- sapply(1:nrow(mat_nn_to), function(j){
  matched_diff <- mat_nn_to[j,] - mat_avg_y[i,]
  stats::cor(pred_diff, matched_diff, method = rec_options$cor_method)
})

plot(cor_vec, col = zz)

list_nn_to[[which.max(cor_vec)]]






