rm(list=ls())
library(Seurat)
library(multiomeFate)

load("../../out/Writeup8b/Writeup8b_simulation_day10-COCL2.RData")

par(mfrow = c(1,1))
Seurat::DimPlot(all_data,
                group.by = "dataset",
                reduction = "umap")

set.seed(10)
embedding_mat <- all_data[["fasttopic_COCL2"]]@cell.embeddings
embedding_mat <- scale(embedding_mat)
embedding_mat <- pmin(embedding_mat, 2)
coefficient_intercept <- -2
embedding_coefficient_vec <- stats::runif(ncol(embedding_mat), min = -0.5, max = 0.25) 

num_lineages <- 100
lineage_prior <- rep(1/num_lineages, length = num_lineages)

############################
# not much to change after this line

early_idx <- which(all_data$dataset == "day10_COCL2")
lineage_spread <- 1
set.seed(10)
simulation_res <- generate_simulation(
  embedding_mat = embedding_mat[early_idx,],
  bool_add_randomness = TRUE,
  coefficient_intercept = coefficient_intercept,
  embedding_coefficient_vec = embedding_coefficient_vec,
  lineage_spread = lineage_spread,
  lineage_prior = lineage_prior,
  num_lineages = num_lineages,
  verbose = 3
)

# table(10^simulation_res$cell_fate_potential)

future_lineage_spread <- 1
future_idx <- which(all_data$dataset == "week5_COCL2")
simulation_future_res <- generate_simulation_attachFuture(
  coefficient_intercept = simulation_res$coefficient_intercept,
  embedding_coefficient_vec = simulation_res$embedding_coefficient_vec,
  future_cell_embedding_mat = embedding_mat[future_idx,],
  lineage_assignment = simulation_res$lineage_assignment,
  previous_cell_embedding_mat = embedding_mat[early_idx,],
  lineage_spread = future_lineage_spread,
  verbose = 2
)


par(mfrow = c(1,2))
plot(all_data[["umap"]]@cell.embeddings,
     pch = 16, col = "gray",
     xlab = "UMAP1", ylab = "UMAP2",
     main = paste0("Spread: ", lineage_spread))
for(kk in 1:3){
  idx <- which(simulation_res$lineage_assignment == paste0("lineage:", kk))
  points(all_data[["umap"]]@cell.embeddings[early_idx[idx],,drop = FALSE],
         pch = 16, col = kk, cex = 1)
}

plot(all_data[["umap"]]@cell.embeddings,
     pch = 16, col = "gray",
     xlab = "UMAP1", ylab = "UMAP2",
     main = paste0("Spread: ", lineage_spread))
future_idx <- which(all_data$dataset == "week5_COCL2")
for(kk in 1:3){
  idx <- which(simulation_res$lineage_assignment == paste0("lineage:", kk))
  previous_cells <- rownames(all_data[["umap"]]@cell.embeddings[early_idx[idx],,drop = FALSE])
  idx2 <- which(simulation_future_res$future_cell_assignment %in% previous_cells)
  
  points(all_data[["umap"]]@cell.embeddings[future_idx[idx2],,drop = FALSE],
         pch = 16, col = "white", cex = 1.5)
  points(all_data[["umap"]]@cell.embeddings[future_idx[idx2],,drop = FALSE],
         pch = 16, col = kk, cex = 1)
}

