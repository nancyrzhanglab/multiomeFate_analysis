rm(list=ls())
library(Seurat)

treatment <- "CIS"
day_early <- "day0"
day_later <- "day10"
day_early_full <- day_early
day_later_full <- paste0(day_later, "_", treatment)

load(paste0("../../../../out/kevin/Writeup6q/Writeup6q_", treatment, "_", day_early, "_lineage-imputation.RData"))

###################

train_mat <- sapply(loocv_fit_list, function(x){
  x$train_loglik
})
Matrix::rowMeans(train_mat)

test_mat <- sapply(loocv_fit_list, function(x){
  x$test_loglik
})
loglik_mean <- Matrix::rowMeans(test_mat)
loglik_mean

lambda_sequence <- loocv_fit_list[[1]]$train_fit$lambda_sequence
lambda <- lambda_sequence[which.min(loglik_mean)]

final_fit <- multiomeFate:::lineage_imputation(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  coefficient_initial_list = loocv_fit_list[[1]]$train_fit$fit_list[[which.min(loglik_mean)]]$coefficient_vec,
  lambda = lambda,
  lineage_future_count = lineage_future_count,
  verbose = 1
)

#######

load("../../../../out/kevin/Writeup6p/Writeup6p_all-data_lightweight_noATAC.RData")

keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = paste0("fasttopic_", treatment),
                            dims = 1:30)

plot1 <- Seurat::DimPlot(all_data, reduction = "umap",
                         group.by = "dataset", pt.size = .3, label = T)
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6q/Writeup6q_", treatment, "_umap.png"),
                plot1, device = "png", width = 5, height = 4, units = "in")

############################

cell_imputed_score <- as.numeric(cell_features %*% final_fit$fit$coefficient_vec)
names(cell_imputed_score) <- rownames(cell_features)
cell_imputed_count <- exp(cell_imputed_score)
uniq_lineage <- sort(unique(cell_lineage))
lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})

lineage_future_count <- tab_mat[uniq_lineage,day_later_full]

###########################

cell_imputed_score_full <- rep(NA, ncol(all_data))
names(cell_imputed_score_full) <- colnames(all_data)

all_data$imputed_count <- cell_imputed_score[colnames(all_data)]
  
p1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        na_cutoff = quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle(paste0(
  treatment, "\n", day_later, " growth potential of ", day_early, 
  " cells\n(Stepdown from RNA fasttopics, ATAC PeakVI)\n(Log-scale)")
)
ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6q/Writeup6q_",
                                  treatment, "-", day_early, "_imputation-ridge_umap.png"),
                p1, device = "png", width = 5, height = 5, units = "in")

###########################

all(names(lineage_imputed_count) == names(lineage_future_count))

lineage_imputed_count2 <- log10(lineage_imputed_count+1)
lineage_future_count2 <- log10(lineage_future_count+1)

labeling_vec <- rep(FALSE, length(lineage_imputed_count2))
labeling_vec[intersect(which(lineage_imputed_count2 >= 1.5),
                       order(lineage_imputed_count2, decreasing = T)[1:10])] <- TRUE
labeling_vec[intersect(which(lineage_future_count2 >= 1.5),
                       order(lineage_future_count2, decreasing = T)[1:10])] <- TRUE

df <- data.frame(lineage_imputed_count = lineage_imputed_count2,
                 lineage_future_count = lineage_future_count2,
                 name = names(lineage_imputed_count2),
                 labeling = labeling_vec)
# put all the labeling == TRUE on bottom
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage_future_count, y = lineage_imputed_count))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = name, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + ggplot2::ggtitle(paste0(
  treatment, " ", day_later, " growth potential of ", day_early, 
  " cells\n(Ridge for RNA fasttopics, ATAC PeakVI), (Log-scale)",
  "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2))
) +
  ggplot2::xlab("Observed lineage count (Log10)") + ggplot2::ylab("Predicted lineage count (Log10)")
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6q/Writeup6q_", treatment, "-", day_early, "_imputation-ridge_lineage-level.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

