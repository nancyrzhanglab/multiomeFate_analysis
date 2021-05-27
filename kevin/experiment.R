rm(list=ls())
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
                           x_sd_biological = 0.5, x_sd_technical = 0.1, 
                           y_exp_baseline = 0.1, y_sd_technical = 1.5,
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
res <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell)
mat_x <- res$obs_x; mat_y <- res$obs_y
idx <- which(df_cell$branch == 1)
mat_x <- mat_x[idx,]; mat_y <- mat_y[idx,]; df_cell <- df_cell[idx,]

image(.rotate(mat_y))
image(.rotate(mat_x))

#############################

idx <- which(df_cell$branch == 1)
mat_x <- mat_x[idx,]; mat_y <- mat_y[idx,]; df_cell <- df_cell[idx,]

vec_start <- which(df_cell$time <= 10)
list_end <- lapply(sort(unique(df_cell$branch)), function(branch){
        intersect(which(df_cell$branch == branch), which(df_cell$time >= 90))
})
# zz <- svd(mat_x); plot(zz$d[1:50])
rank_x <- 20
# zz <- svd(mat_y); plot(zz$d[1:50])
rank_y <- 20
set.seed(10)
prep_obj <- chromatin_potential_prepare(mat_x, mat_y, df_x, df_y, 
                                        vec_start, list_end,
                                        est_method = "threshold_glmnet",
                                        options = list(nn_nn = 10, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, est_cis_window = 30,
                                                       est_num_iterations = 4))

set.seed(10)
verbose = T
mat_g_init = NA; vec_g_init = rep(0, ncol(mat_y))

#####################################
mat_x <- prep_obj$mat_x; mat_y <- prep_obj$mat_y
df_x <- prep_obj$df_x; df_y <- prep_obj$df_y
df_res <- prep_obj$df_res; dim_reduc_obj <- prep_obj$dim_reduc_obj
nn_mat <- prep_obj$nn_mat; nn_obj <- prep_obj$nn_obj
list_diagnos <- prep_obj$list_diagnos; options <- prep_obj$options

dim_options <- options$dim_options; nn_options <- options$nn_options
form_options <- options$form_options; est_options <- options$est_options
cand_options <- options$cand_options; rec_options <- options$rec_options

# initialize
n <- nrow(mat_x)

ht_neighbor <- .init_chrom_ht(which(df_res$order_rec == 0))
tmp <- .init_est_matrices(mat_x, mat_y, df_res)
mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
list_diagnos <- list()
iter <- 1

if(verbose) print(paste0("Iteration ", iter, ": Recruited percentage (", 
                         round(sum(!is.na(df_res$order_rec))/nrow(df_res), 2), ")"))
## estimate res_g
if((iter == 1 | est_options$hold_initial) && !any(is.na(mat_g_init)) && !any(is.na(vec_g_init))){
        res_g <- list(mat_g = mat_g_init, vec_g = vec_g_init)
} else {
        res_g <- .estimate_g(mat_x1, mat_y2, est_options)
}

## construct candidate set
res_cand <- .candidate_set(mat_x, mat_y, df_res, nn_mat, cand_options)
df_res <- .update_chrom_df_cand(df_res, res_cand$vec_cand)
stopifnot(all(is.na(df_res$order_rec[res_cand$vec_cand])))
list_diagnos[[as.character(iter)]]$candidate <- res_cand$diagnostic

## recruit an element from the candidate set
enforce_matched <- length(which(df_res$order_rec == 0)) > length(which(df_res$order_rec > 0))
res_rec <- .recruit_next(mat_x, mat_y, res_cand$vec_cand, res_g, df_res, 
                         dim_reduc_obj, nn_mat, nn_obj, enforce_matched,
                         rec_options)
stopifnot(all(is.na(df_res$order_rec[res_rec$rec$vec_from])))
list_diagnos[[as.character(iter)]]$recruit <- res_rec$diagnostic

## update
tmp <- .update_estimation_matrices(mat_x, mat_y, mat_x1, mat_y2, 
                                   res_rec$rec, form_options)
mat_x1 <- tmp$mat_x1; mat_y2 <- tmp$mat_y2
ht_neighbor <- .update_chrom_ht(ht_neighbor, res_rec$rec$vec_from, 
                                res_rec$rec$list_to, enforce_matched)
df_res <- .update_chrom_df_rec(df_res, res_rec$rec$vec_from, iter)

iter <- iter+1

####

if(verbose) print(paste0("Iteration ", iter, ": Recruited percentage (", 
                         round(sum(!is.na(df_res$order_rec))/nrow(df_res), 2), ")"))
## estimate res_g
if((iter == 1 | est_options$hold_initial) && !any(is.na(mat_g_init)) && !any(is.na(vec_g_init))){
        res_g <- list(mat_g = mat_g_init, vec_g = vec_g_init)
} else {
        res_g <- .estimate_g(mat_x1, mat_y2, est_options)
}

## construct candidate set
res_cand <- .candidate_set(mat_x, mat_y, df_res, nn_mat, cand_options)
df_res <- .update_chrom_df_cand(df_res, res_cand$vec_cand)
stopifnot(all(is.na(df_res$order_rec[res_cand$vec_cand])))
list_diagnos[[as.character(iter)]]$candidate <- res_cand$diagnostic

########
image(.rotate(abs(res_g$mat_g)))
plot(res_g$vec_g)
plot(res_g$vec_threshold)

########
vec_cand <- res_cand$vec_cand
pred_y <- .predict_yfromx(mat_x[vec_cand,,drop = F], res_g, rec_options$family)

nn_size <- ncol(nn_mat)
i <- 1
cell <- vec_cand[i]

nn_cand <- c(nn_mat[cell, ], cell)

vec <- c(.apply_dimred(mat_x[vec_cand[i],], dim_reduc_obj$x),
         .apply_dimred(pred_y[i,], dim_reduc_obj$y))

nn_pred <- .distant_nn(cell, nn_mat)

# find all nn's that aren't too close to cell itself
if(length(setdiff(nn_pred, nn_cand)) > 0) nn_pred <- setdiff(nn_pred, nn_cand)
df_cell[nn_pred,]

# from this set of cells, find the ones with highest pearson
# [[note to self: this should be refactored out]]
pred_diff <- pred_y[i,] - mat_y[vec_cand[i],]
cor_vec <- sapply(nn_pred, function(j){
        matched_diff <- mat_y[j,] - mat_y[vec_cand[i],]
        stats::cor(pred_diff, matched_diff, method = rec_options$cor_method)
})
idx <- nn_pred[which.max(cor_vec)]

#######



xlim <- c(1, ncol(mat_y))
ylim <- range(c(pred_y, mat_y))
plot(NA, xlim = xlim, ylim = ylim)
for(i in 1:nrow(pred_y)){
        points(pred_y[i,], pch = 16)
}
for(i in res_cand$vec_cand){
        points(mat_y[i,], pch = 16, col = "red")
}

plot(pred_y[1,], mat_y[1,], asp = T)
