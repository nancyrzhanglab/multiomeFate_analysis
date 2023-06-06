rm(list=ls())

set.seed(10)
n_each <- 30
res <- .construct_lineage_data(coefficient_vec = c(1,0),
                               L = 10, 
                               n_each = n_each, 
                               variance_within_lineage = 0.5)
cell_features <- res$cell_features
cell_lineage <- res$cell_lineage
cell_lineage_idx_list <- res$cell_lineage_idx_list
true_coefficient <- res$coefficient_vec
coefficient_initial <- true_coefficient/2
lineage_future_count <- res$lineage_future_count

res <- lineage_imputation(cell_features,
                          cell_lineage,
                          coefficient_initial,
                          lineage_future_count)

other_obj <- .lineage_objective(cell_features = cell_features,
                                cell_lineage = cell_lineage,
                                cell_lineage_idx_list = cell_lineage_idx_list,
                                coefficient_vec = true_coefficient,
                                lineage_future_count = lineage_future_count)

# make two plots, one of raw data, the other of the imputed data
par(mfrow = c(1,3), mar = rep(0.5,4))
uniq_lineage <- sort(unique(cell_lineage))
col_palette <- scales::hue_pal()(length(uniq_lineage))
names(col_palette) <- uniq_lineage
col_vec <- col_palette[cell_lineage]
size_vec <- lineage_future_count[cell_lineage]/(50*n_each)
plot(cell_features[,1], cell_features[,2], col = col_vec, cex = size_vec, pch = 16)

imputed_counts <- as.numeric(exp(cell_features %*% res$coefficient_vec))
names(imputed_counts) <- rownames(cell_features)
size_vec2 <- imputed_counts/25
size_vec2 <- pmin(size_vec2, quantile(size_vec2, prob = 0.99))
plot(cell_features[,1], cell_features[,2],
     col = col_vec, cex = size_vec2, pch = 16)


true_counts <- as.numeric(exp(cell_features %*% true_coefficient))
names(true_counts) <- rownames(cell_features)
size_vec3 <- true_counts/25
size_vec3 <- pmin(size_vec3, quantile(size_vec3, prob = 0.99))
plot(cell_features[,1], cell_features[,2],
     col = col_vec, cex = size_vec3, pch = 16)

par(mfrow = c(1,2), mar = c(4,4,4,0.5))
plot(true_counts, imputed_counts, pch = 16, asp = T)

lineage_true_counts <- sapply(cell_lineage_idx_list, function(idx){
  sum(true_counts[idx])
})
lineage_imputed_counts <- sapply(cell_lineage_idx_list, function(idx){
  sum(imputed_counts[idx])
})

plot(lineage_true_counts, lineage_imputed_counts, pch = 16, asp = T)


