rm(list=ls())
trial <- 1
res <- .construct_lineage_data(p = 1, seed = trial)
cell_features <- res$cell_features
cell_lineage <- res$cell_lineage
cell_lineage_idx_list <- res$cell_lineage_idx_list
lineage_future_count <- res$lineage_future_count

set.seed(trial)
coef_target <- runif(1)

obj_val <- .lineage_objective(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              cell_lineage_idx_list = cell_lineage_idx_list,
                              coefficient_vec = coef_target,
                              lineage_future_count = lineage_future_count)
grad_val <- .lineage_gradient(cell_features = cell_features,
                              cell_lineage = cell_lineage,
                              cell_lineage_idx_list = cell_lineage_idx_list,
                              coefficient_vec = coef_target,
                              lineage_future_count = lineage_future_count)

coef_jitter <- coef_target + seq(-10,2.25,by=0.01)
obj_vec <- sapply(coef_jitter, function(coef_val){
  .lineage_objective(cell_features = cell_features,
                     cell_lineage = cell_lineage,
                     cell_lineage_idx_list = cell_lineage_idx_list,
                     coefficient_vec = coef_val,
                     lineage_future_count = lineage_future_count)
})

lower_bound_vec <- sapply(coef_jitter, function(coef_val){
  obj_val - grad_val*(coef_target-coef_val)
})

plot(coef_jitter, obj_vec)
points(coef_target, obj_val, pch = 16, col = "red")
lines(coef_jitter, lower_bound_vec, col = "red")
