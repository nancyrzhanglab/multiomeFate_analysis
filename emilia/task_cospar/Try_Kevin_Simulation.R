library(multiomeFate)
library(Seurat)

load("~/Downloads/Writeup8b_simulation_day10-COCL2.RData")

par(mfrow = c(1,1))
Seurat::DimPlot(all_data,
                group.by = "dataset",
                reduction = "umap")

embedding_mat <- all_data[["fasttopic_COCL2"]]@cell.embeddings
embedding_mat <- scale(embedding_mat)
embedding_mat <- pmin(embedding_mat, 2)
coefficient_intercept <- -1
coefficient_vec <- rep(1, ncol(embedding_mat))
coefficient_vec <- coefficient_vec/2

# double-check the fate potentials are ok
tmp <- exp((embedding_mat %*% coefficient_vec) + coefficient_intercept)
quantile(tmp)
sum(tmp)

lineage_spread <- 100
num_lineages <- 50
lineage_prior <- rep(1/num_lineages, length = num_lineages)

############################
# not much to change after this line

early_idx <- which(all_data$dataset == "day10_COCL2")
set.seed(10)
simulation_res <- multiomeFate:::generate_simulation(
  embedding_mat = embedding_mat[early_idx,,drop=FALSE],
  coefficient_intercept = coefficient_intercept,
  # coefficient_vec = coefficient_vec,
  lineage_spread = lineage_spread,
  lineage_prior = lineage_prior,
  num_lineages = num_lineages,
  verbose = 3
)

# check the simulation to make the sizes look alright
table(simulation_res$lineage_assignment)
hist(simulation_res$cell_fate_potential)
hist(10^(simulation_res$cell_fate_potential))
hist(simulation_res$lineage_future_size)
sum(simulation_res$lineage_future_size)
sum(10^(simulation_res$cell_fate_potential))
round(simulation_res$prob_mat*num_lineages, 2)

par(mfrow = c(1,1))
idx <- which(simulation_res$lineage_assignment == "lineage:1")
plot(all_data[["umap"]]@cell.embeddings,
     pch = 16, col = "gray",
     xlab = "UMAP1", ylab = "UMAP2",
     main = paste0("Spread: ", lineage_spread))
for(i in idx){
  points(all_data[["umap"]]@cell.embeddings[early_idx[idx],,drop = FALSE],
         pch = 16, col = "coral3", cex = 1)
}

##################

# try out our fate potential method
cell_lineage <- as.character(simulation_res$lineage_assignment)
uniq_lineage <- sort(unique(cell_lineage))
lineage_future_count <- simulation_res$lineage_future_size
tmp <- table(simulation_res$lineage_assignment)
lineage_current_count <- as.numeric(tmp); names(lineage_current_count) <- names(tmp)
lineage_current_count <- lineage_current_count[names(lineage_future_count)]
tab_mat <- cbind(lineage_current_count, lineage_future_count)
colnames(tab_mat) <- c("now", "future")

#################
# start cross validation

fit_res <- multiomeFate:::lineage_cv(
  cell_features = embedding_mat[early_idx,,drop=FALSE],
  cell_lineage = cell_lineage,
  future_timepoint = "future",
  lineage_future_count = lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 10,
  tab_mat = tab_mat,
  num_folds = 10,
  verbose = 2
)

final_fit <- multiomeFate:::lineage_cv_finalize(
  cell_features = embedding_mat[early_idx,,drop=FALSE],
  cell_lineage = cell_lineage,
  fit_res = fit_res,
  lineage_future_count = lineage_future_count
)
lineage_imputed_count <- final_fit$lineage_imputed_count
cell_imputed_score <- final_fit$cell_imputed_score
round(final_fit$coefficient_vec, 2)

par(mfrow = c(1,2))
plot(x = simulation_res$lineage_future_size,
     y = lineage_imputed_count[names(simulation_res$lineage_future_size)], 
     asp = T, pch = 16,
     xlab = "True lineage size",
     ylab = "Fitted lineage size",
     main = paste0("Correlation: ", round(stats::cor(
       simulation_res$lineage_future_size,
       lineage_imputed_count[names(simulation_res$lineage_future_size)]
     ), 2))
)

plot(x = simulation_res$cell_fate_potential,
     y = cell_imputed_score[names(simulation_res$cell_fate_potential)], 
     asp = T, pch = 16,
     xlab = "True cell potential",
     ylab = "Estimated cell potential",
     main = paste0("Correlation: ", round(stats::cor(
       simulation_res$cell_fate_potential,
       cell_imputed_score[names(simulation_res$cell_fate_potential)]
     ), 2))
)
