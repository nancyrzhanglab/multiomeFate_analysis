# this simulation is for similar means but the different ranges

rm(list=ls())
library(Seurat)
library(multiomeFate)

load("~/project/Multiome_fate/out/kevin/Writeup9b/Writeup9b_simulation_day10-COCL2.RData")
load("~/project/Multiome_fate/out/kevin/Writeup9b/Writeup9b_day10-COCL2_saver.RData")

set.seed(10)
rna_mat <- saver_res$estimate
rna_mat <- Matrix::t(rna_mat)
rna_mat <- pmin(rna_mat, 10)

early_idx <- which(all_data$dataset == "day10_COCL2")
rna_mat <- rna_mat[early_idx,]

mean_vec <- colMeans(rna_mat)
sd_vec <- apply(rna_mat, 2, stats::sd)
rm_idx <- which(sd_vec <= 1e-3)
if(length(rm_idx) > 0){
  rna_mat <- rna_mat[,-rm_idx,drop = FALSE]
}
cell_features <- scale(rna_mat)
svd_res <- irlba::irlba(cell_features, nv = 30)
svd_features <- sweep(svd_res$u, MARGIN = 2, STATS = svd_res$d, FUN = '*')
rownames(svd_features) <- rownames(rna_mat)
colnames(svd_features) <- paste0("svd_", 1: ncol(svd_features))
cell_features <- svd_features

coefficient_intercept <- -14
svd_coefficient_vec <- rep(0, ncol(svd_features))
names(svd_coefficient_vec) <- colnames(svd_features)
svd_coefficient_vec[1:5] <- seq(0.2, 0.1, length.out = 5)
coefficient_vec <- svd_coefficient_vec
names(coefficient_vec) <- colnames(cell_features)

# double-check the fate potentials are ok
tmp2 <- exp((cell_features %*% coefficient_vec) + coefficient_intercept)
sum(tmp2)
# plot(tmp, tmp2, main = paste0("Corr: ", round(cor(tmp, tmp2), 2)))

num_lineages <- 50

############################
# not much to change after this line
set.seed(10)
simulation_res <- multiomeFate:::generate_simulation_plastic(
  embedding_mat = cell_features,
  bool_add_randomness = TRUE,
  coefficient_intercept = coefficient_intercept,
  embedding_coefficient_vec = coefficient_vec,
  num_lineages = num_lineages,
  lineage_mean_spread = 1, 
  lineage_sd_spread = NA,
  verbose = 3
)

# check the simulation to make the sizes look alright
# table(simulation_res$lineage_assignment)
# hist(simulation_res$cell_fate_potential)
# hist(10^(simulation_res$cell_fate_potential)-1)
# hist(simulation_res$lineage_future_size)
# sum(simulation_res$lineage_future_size)
# sum(10^(simulation_res$cell_fate_potential))
# 
# par(mfrow = c(1,1))
# plot(simulation_res$cell_fate_potential_truth,
#      simulation_res$cell_fate_potential,
#      pch = 16, asp = TRUE)
# 
# par(mfrow = c(1,1))
# plot(all_data[["umap"]]@cell.embeddings,
#      pch = 16, col = "gray",
#      xlab = "UMAP1", ylab = "UMAP2")
# lineage_idx_vec <- c(1,49,50)
# for(kk in 1:length(lineage_idx_vec)){
#   idx <- which(simulation_res$lineage_assignment == paste0("lineage:", lineage_idx_vec[kk]))
#   points(all_data[["umap"]]@cell.embeddings[early_idx[idx],,drop = FALSE],
#          pch = 16, col = kk, cex = 1)
# }
# 
# par(mfrow = c(1,1))
# plot(x = simulation_res$cell_fate_potential_truth,
#      y = simulation_res$lineage_future_size[simulation_res$lineage_assignment],
#      xlab = "Cell GP", ylab = "Lineage future size",
#      pch = 16, cex = 0.5)

png(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9b/Writeup9b_alldata-cocl2_simulation_variance2growth.png"),
    width = 15, height = 5, units = "in", res = 300)
par(mfrow = c(1,3))
median_vec <- simulation_res$summary_mat["median",]
range_vec <- simulation_res$summary_mat["range",]
future_vec <- simulation_res$summary_mat["future_size",]
plot(median_vec, 
     future_vec, 
     pch = 16,
     xlab = "Median (per lineage)",
     ylab = "Future size (per lineage)",
     main = paste0("Future size vs. current median potential\n",
                   "Corr: ", round(cor(median_vec, future_vec), 2)))
plot(range_vec,
     future_vec,
     pch = 16,
     xlab = "Range (per lineage)",
     ylab = "Future size (per lineage)",
     main = paste0("Future size vs. current range potential\n",
                   "Corr: ", round(cor(range_vec, future_vec), 2)))
plot(x = median_vec,
     y = range_vec, 
     main = paste0("Corr: ", round(cor(range_vec, median_vec), 2)),
     xlab = "Median (per lineage)", ylab = "Range (per lineage)",
     pch = 16, asp = TRUE)
graphics.off()

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

set.seed(10)
start_time1 <- Sys.time()
fit_res <- multiomeFate:::lineage_cv(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  future_timepoint = "future",
  lineage_future_count = lineage_future_count,
  lambda_initial = 3,
  lambda_sequence_length = 10,
  tab_mat = tab_mat,
  num_folds = 10,
  verbose = 2
)
end_time1 <- Sys.time()

start_time2 <- Sys.time()
final_fit <- multiomeFate:::lineage_cv_finalize(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  fit_res = fit_res,
  lineage_future_count = lineage_future_count
)
end_time2 <- Sys.time()
lineage_imputed_count <- final_fit$lineage_imputed_count
cell_imputed_score <- final_fit$cell_imputed_score
round(final_fit$coefficient_vec, 2)

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(fit_res, final_fit, simulation_res,
     date_of_run, session_info,
     start_time1, start_time2, end_time1, end_time2,
     file = "~/project/Multiome_fate/out/kevin/Writeup9b/Writeup9b_day10-COCL2_saver-pca_simulation_variance2growth.RData")

print("Done! :)")
