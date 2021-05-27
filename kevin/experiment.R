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
res <- simulate_data(df_x, df_y, list_xnoise, list_ynoise, df_cell,
                     jitter = 0.1)
mat_x <- res$mean_x; mat_y <- res$mean_y

###############################

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
                                        rec_method = "distant_cor_oracle",
                                        options = list(nn_nn = 10, dim_nlatent_x = rank_x,
                                                       dim_nlatent_y = rank_y, est_cis_window = 30,
                                                       est_num_iterations = 4))

set.seed(10)
# res <- chromatin_potential(prep_obj, df_cell = df_cell, verbose = T)

################################
mat_g_init = NA
vec_g_init = rep(0, ncol(mat_y))
verbose = T

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
# res_g <- .estimate_g(mat_x1, mat_y2, est_options)

########################################

stopifnot(nrow(mat_x1) == nrow(mat_y2))
if(est_options$enforce_cis) stopifnot(class(est_options$ht_map) == "hash")

p1 <- ncol(mat_x1); p2 <- ncol(mat_y2)

# initialize variables for the loop
if(!est_options$parallel && future::nbrOfWorkers() == 1){
        my_lapply <- pbapply::pblapply
        if(est_options$verbose) pbapply::pboptions(type = "timer") else pbapply::pboptions(type = "none")
} else {
        my_lapply <- future.apply::future_lapply
}

j <- 118
set.seed(10)
print(j)
if(est_options$enforce_cis){
        ## find the region around each peak
        idx_x <- est_options$ht_map[[as.character(j)]]
} else {
        #[[note to self: I think there's a cleaner way to write this]]
        #[[note to self: write a test for this scenario]]
        idx_x <- 1:p1
}
if(length(idx_x) == 0) return(val_int = mean(mat_y2[,j]), vec_coef = rep(0, p1))

idx_y <- j

print(quantile(mat_y2[,idx_y]))

## apply glmnet
res <- .threshold_glmnet_fancy(mat_x1[,idx_x,drop = F], mat_y2[,idx_y],
                        family = est_options$family, 
                        switch = est_options$switch, switch_cutoff = est_options$switch_cutoff,
                        alpha = est_options$alpha, standardize = est_options$standardize, intercept = est_options$intercept,
                        cv = est_options$cv, nfolds = est_options$nfolds, cv_choice = est_options$cv_choice,
                        bool_round = est_options$bool_round,
                        num_iterations = est_options$num_iterations, 
                        initial_quantile = est_options$initial_quantile)


#.transform_est_matrix(list_res, est_options, p1)
