rm(list=ls())

set.seed(10)
res <- .construct_lineage_data()
cell_features <- res$cell_features
cell_lineage <- res$cell_lineage
cell_lineage_idx_list <- res$cell_lineage_idx_list
true_coefficient <- res$coefficient_vec
coefficient_initial <- true_coefficient/2
lineage_future_count <- res$lineage_future_count

lineage_imputation(cell_features,
                   cell_lineage,
                   coefficient_initial,
                   lineage_future_count)

.lineage_objective(cell_features = cell_features,
                   cell_lineage = cell_lineage,
                   cell_lineage_idx_list = cell_lineage_idx_list,
                   coefficient_vec = true_coefficient,
                   lineage_future_count = lineage_future_count)
