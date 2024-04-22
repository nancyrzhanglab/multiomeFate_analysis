rm(list=ls())
library(Seurat)
library(multiomeFate)

load("../../out/Writeup8b/Writeup8b_simulation_day10-COCL2.RData")

# par(mfrow = c(1,1))
# Seurat::DimPlot(all_data,
#                 group.by = "dataset",
#                 reduction = "umap")

set.seed(10)
embedding_mat <- all_data[["fasttopic_COCL2"]]@cell.embeddings
embedding_mat <- scale(embedding_mat)
coefficient_intercept <- -8
embedding_coefficient_vec <- rep(0, ncol(embedding_mat))
names(embedding_coefficient_vec) <- colnames(embedding_mat)
embedding_coefficient_vec["fastTopicCOCL2_1"] <- 2

held_out_variables <- c("fastTopicCOCL2_2", "fastTopicCOCL2_3", "fastTopicCOCL2_4")
embedding_coefficient_vec[held_out_variables] <- 1.1 # this is on the exponential scale, so it's **very** sensitive

rm_idx <- which(colnames(embedding_mat) %in% held_out_variables)
fatefeatures_mat <- embedding_mat[,rm_idx,drop = FALSE]
fatefeatures_coefficient_vec <- embedding_coefficient_vec[rm_idx]

embedding_mat <- embedding_mat[,-rm_idx,drop = FALSE]
embedding_coefficient_vec <- embedding_coefficient_vec[-rm_idx]

# double-check the fate potentials are ok
tmp <- exp((embedding_mat %*% embedding_coefficient_vec) + coefficient_intercept)
quantile(tmp)
tmp2 <- exp((embedding_mat %*% embedding_coefficient_vec) + (fatefeatures_mat %*% fatefeatures_coefficient_vec) + coefficient_intercept)
sum(tmp2)
cor(tmp, tmp2) # the lower this number is, the harder the problem is
plot(tmp, tmp2)

num_lineages <- 50

############################
# not much to change after this line

early_idx <- which(all_data$dataset == "day10_COCL2")
set.seed(10)
simulation_res <- multiomeFate:::generate_simulation_plastic(
  embedding_mat = embedding_mat[early_idx,,drop=FALSE],
  bool_add_randomness = TRUE,
  coefficient_intercept = coefficient_intercept,
  embedding_coefficient_vec = embedding_coefficient_vec,
  fatefeatures_coefficient_vec = fatefeatures_coefficient_vec,
  fatefeatures_mat = fatefeatures_mat[early_idx,,drop=FALSE], 
  num_lineages = num_lineages,
  verbose = 3
)

# check the simulation to make the sizes look alright
table(simulation_res$lineage_assignment)
hist(simulation_res$cell_fate_potential)
hist(10^(simulation_res$cell_fate_potential)-1)
hist(simulation_res$lineage_future_size)
sum(simulation_res$lineage_future_size)
sum(10^(simulation_res$cell_fate_potential))

par(mfrow = c(1,1))
plot(simulation_res$cell_fate_potential_truth,
     simulation_res$cell_fate_potential,
     pch = 16, asp = TRUE)

par(mfrow = c(1,1))
plot(all_data[["umap"]]@cell.embeddings,
     pch = 16, col = "gray",
     xlab = "UMAP1", ylab = "UMAP2")
lineage_idx_vec <- c(1,49,50)
for(kk in 1:length(lineage_idx_vec)){
  idx <- which(simulation_res$lineage_assignment == paste0("lineage:", lineage_idx_vec[kk]))
  points(all_data[["umap"]]@cell.embeddings[early_idx[idx],,drop = FALSE],
         pch = 16, col = kk, cex = 1)
}

par(mfrow = c(1,3))
mean_val = sapply(levels(simulation_res$lineage_assignment), function(lineage){
  median(simulation_res$cell_fate_potential_truth[simulation_res$lineage_assignment == lineage])
})
plot(mean_val, 
     simulation_res$lineage_future_size, 
     pch = 16,
     xlab = "Mean (per lineage)",
     ylab = "Future size (per lineage)",
     main = paste0("Future size vs. current mean potential\n",
                   "Corr: ", round(cor(mean_val, simulation_res$lineage_future_size), 2)))
sd_val = sapply(levels(simulation_res$lineage_assignment), function(lineage){
  diff(range(simulation_res$cell_fate_potential_truth[simulation_res$lineage_assignment == lineage]))
})
plot(sd_val, 
     simulation_res$lineage_future_size, 
     pch = 16,
     xlab = "Sd (per lineage)",
     ylab = "Future size (per lineage)",
     main = paste0("Future size vs. current sd potential\n",
                   "Corr: ", round(cor(sd_val, simulation_res$lineage_future_size), 2)))
plot(x = mean_val,
     y = sd_val, 
     main = paste0("Corr: ", round(cor(sd_val, mean_val), 2)),
     xlab = "Mean (per lineage)", ylab = "Sd (per lineage)",
     pch = 16, asp = TRUE)

#####

lineage_vec <- simulation_res$lineage_assignment
later_size <- simulation_res$lineage_future_size
lineage_names <- names(sort(later_size, decreasing = T))[1:20]
idx <- which(lineage_vec %in% lineage_names)

df <- data.frame(lineage = lineage_vec[idx],
                 imputed_count = simulation_res$cell_fate_potential_truth[idx])
df_tmp <- df; df_tmp$lineage <- droplevels(as.factor(df_tmp$lineage))
anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df_tmp)
df2 <- data.frame(lineage = "All",
                  imputed_count = simulation_res$cell_fate_potential_truth)
df <- rbind(df, df2)

total_std <- sum((df_tmp$imputed_count - mean(df_tmp$imputed_count))^2)
within_lineage_std <- sum(sapply(levels(df_tmp$lineage), function(lineage_name){
  idx <- which(df_tmp$lineage == lineage_name)
  sum((df_tmp$imputed_count[idx] - mean(df_tmp$imputed_count[idx]))^2)
}))
across_lineage_std <- sum(sapply(levels(df_tmp$lineage), function(lineage_name){
  idx <- which(df_tmp$lineage == lineage_name)
  mean_val <- mean(df_tmp$imputed_count[idx])
  length(idx) * (mean_val - mean(df_tmp$imputed_count))^2 
}))
lineage_effect <- round(across_lineage_std/total_std*100,1)

col_vec <- c(rep("#999999", length(lineage_names)), "#E69F00")
names(col_vec) <- c(lineage_names, "All")

p1 <- ggplot2::ggplot(df, ggplot2::aes(x=lineage, y=imputed_count))
p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width", ggplot2::aes(fill=lineage))
p1 <- p1 + ggplot2::scale_fill_manual(values = col_vec) 
p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
p1 <- p1 + ggplot2::scale_x_discrete(limits = c(lineage_names, "All"),
                                     guide = ggplot2::guide_axis(angle = 45))
p1 <- p1 + ggplot2::ylab("Week5 growth potential")
p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
p1 <- p1 + ggplot2::stat_summary(fun=max, geom="point", shape=10, size=5, color="blue")
p1 <- p1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2), ", Lineage effect = ", lineage_effect, "%"))
p1

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

par(mfrow = c(1,1))
plot(x = c(coefficient_intercept, embedding_coefficient_vec),
     y = final_fit$coefficient_vec,
     xlab = "True coefficients",
     ylab = "Estimated coefficients",
     asp = TRUE, pch = 16)

par(mfrow = c(1,2))
plot(x = simulation_res$lineage_future_size,
     y = lineage_imputed_count[names(simulation_res$lineage_future_size)], 
     asp = TRUE, pch = 16,
     xlab = "Observed lineage size",
     ylab = "Fitted lineage size",
     main = paste0("Correlation: ", round(stats::cor(
       simulation_res$lineage_future_size,
       lineage_imputed_count[names(simulation_res$lineage_future_size)]
     ), 2))
)

plot(x = simulation_res$cell_fate_potential_truth,
     y = cell_imputed_score[names(simulation_res$cell_fate_potential_truth)], 
     asp = TRUE, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.3),
     xlab = "True cell potential",
     ylab = "Estimated cell potential",
     main = paste0("Correlation: ", round(stats::cor(
       simulation_res$cell_fate_potential_truth,
       cell_imputed_score[names(simulation_res$cell_fate_potential_truth)]
     ), 2))
)



