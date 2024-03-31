rm(list=ls())
library(Seurat)
library(multiomeFate)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step3_fasttopics.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

seurat_object_safe <- seurat_object

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  
  seurat_object <- seurat_object_safe
  keep_vec <- rep(FALSE, ncol(seurat_object))
  idx <- which(seurat_object$time_celltype %in% c(day_early_vec,treatment))
  keep_vec[idx] <- TRUE
  seurat_object$keep <- keep_vec
  seurat_object <- subset(seurat_object, keep == TRUE)
  
  seurat_object <- Seurat::RunUMAP(seurat_object,
                                   dims = 1:30,
                                   reduction = "fasttopic",
                                   reduction.name = "ft_umap")
  
  day_later <- treatment
  
  plot1 <- Seurat::DimPlot(seurat_object, reduction = "ft_umap",
                           group.by = "time_celltype", pt.size = .3, label = T)
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_", treatment, "_umap.png"),
                  plot1, device = "png", width = 6, height = 4, units = "in")
  
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_", treatment, "_from_day", day_early, "_lineage-imputation.RData"))
  
  ###################
  
  train_mat <- sapply(cv_fit_list, function(x){
    x$train_loglik
  })
  
  test_mat <- sapply(cv_fit_list, function(x){
    x$test_loglik
  })
  
  train_quantile <- apply(train_mat, 1, function(vec){stats::quantile(vec, probs = c(0.1, 0.5, 0.9))})
  test_quantile <- apply(test_mat, 1, function(vec){stats::quantile(vec, probs = c(0.1, 0.5, 0.9))})
  
  lambda_sequence <- cv_fit_list[[1]]$train_fit$lambda_sequence
  lambda <- lambda_sequence[which.min(test_quantile[2,])]
  lambda
  
  final_fit <- multiomeFate:::lineage_imputation(
    cell_features = cell_features,
    cell_lineage = cell_lineage,
    coefficient_initial_list = cv_fit_list[[1]]$train_fit$fit_list[[which.min(test_quantile[2,])]]$coefficient_vec,
    lambda = lambda,
    lineage_future_count = lineage_future_count,
    verbose = 0
  )
  
  #######
  
  print("Plotting train-test curve")
  
  lambda_sequence2 <- lambda_sequence+1
  max_base10 <- max(ceiling(log10(lambda_sequence+1)))
  xaxt_values <- unlist(lapply(1:max_base10, function(x){
    seq(10^(x-1), 10^x, length.out = 10)
  }))
  xaxt_values <- xaxt_values[!duplicated(xaxt_values)]
  
  png(paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_",
             treatment, "_train-test-curve.png"),
      height = 1500, width = 3000, units = "px", res = 300)
  par(mfrow = c(1,2), mar = c(4,4,4,4))
  plot(lambda_sequence2, train_quantile[2,], 
       log = "x", type = "n",
       main = paste0(treatment, " ", day_later, " growth potential of ", day_early, "\n(Training)"),
       xlab = "Lambda+1 (Log-scale tickmarks)", 
       ylab = "Negative loglikelihood (Training)",
       ylim = range(c(train_quantile[1,], rev(train_quantile[3,]))))
  polygon(x = c(lambda_sequence2, rev(lambda_sequence2)),
          y = c(train_quantile[1,], rev(train_quantile[3,])),
          col = "gray", 
          border = "black")
  points(lambda_sequence2, train_quantile[2,], pch = 16)
  lines(lambda_sequence2, train_quantile[2,], lwd = 2)
  
  test_upper <- test_quantile[3,]
  test_upper <- pmin(test_upper, stats::quantile(test_upper, probs = 0.95))
  test_lower <- test_quantile[1,]
  test_upper <- pmax(test_upper, stats::quantile(test_upper, probs = 0.05))
  
  plot(lambda_sequence2, test_quantile[2,], 
       log = "x", type = "n",
       main = paste0(treatment, " ", day_later, " growth potential of ", day_early,
                     "\nLambda = ", round(lambda,3), " (Testing)"),
       xlab = "Lambda+1 (Log-scale tickmarks)", 
       ylab = "Negative loglikelihood (Testing)",
       ylim = range(c(test_upper, rev(test_lower))))
  polygon(x = c(lambda_sequence2, rev(lambda_sequence2)),
          y = c(test_upper, rev(test_lower)),
          col = "gray", 
          border = "black")
  points(lambda_sequence2, test_quantile[2,], pch = 16)
  lines(lambda_sequence2, test_quantile[2,], lwd = 2)
  lines(rep(lambda+1, 2), c(-1e6, 1e6), col = "red", lty = 2)
  
  graphics.off()
  
  ############################
  
  cell_features <- cbind(1, cell_features)
  colnames(cell_features)[1] <- "Intercept"
  stopifnot(all(colnames(cell_features) == names(final_fit$fit$coefficient_vec)))
  cell_imputed_score <- as.numeric(cell_features %*% final_fit$fit$coefficient_vec)
  names(cell_imputed_score) <- rownames(cell_features)
  cell_imputed_count <- exp(cell_imputed_score)
  uniq_lineage <- sort(unique(cell_lineage))
  lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
    sum(cell_imputed_count[which(cell_lineage == lineage)])
  })
  cell_imputed_score2 <- log10(exp(cell_imputed_score)) # this one is on the log10 scale
  
  lineage_future_count <- tab_mat[uniq_lineage,day_later]
  
  ###########################
  
  lineage_vec <- seurat_object$assigned_lineage[names(cell_imputed_score2)]
  tab_mat <- table(seurat_object$assigned_lineage, seurat_object$time_celltype)
  tab_vec <- table(lineage_vec)
  tab_vec <- tab_vec[tab_vec >= 2]
  later_size <- tab_mat[names(tab_vec), day_later]
  lineage_names <- names(sort(later_size, decreasing = T))[1:20]
  idx <- which(lineage_vec %in% lineage_names)
  
  df <- data.frame(lineage = lineage_vec[idx],
                   imputed_count = cell_imputed_score2[idx])
  df_tmp <- df; df_tmp$lineage <- as.factor(df_tmp$lineage)
  anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    imputed_count = cell_imputed_score2)
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
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_",
                                    treatment, "_lineage-growthpotential-violinplot.png"),
                  p1, device = "png", width = 6, height = 3, units = "in")
  
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  
  cell_imputed_score_full <- rep(NA, ncol(seurat_object))
  names(cell_imputed_score_full) <- colnames(seurat_object)
  cell_imputed_score_full[names(cell_imputed_score2)] <- cell_imputed_score2
  
  seurat_object$imputed_count <- cell_imputed_score_full
  max_val <- stats::quantile(cell_imputed_score_full, probs = 0.99, na.rm = T)
  seurat_object$imputed_count_thres <- pmin(cell_imputed_score_full, max_val)
  
  p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                          colors_use = list("red", "lightgray", "blue"),
                                          na_cutoff = quantile(seurat_object$imputed_count_thres, probs = 0.05, na.rm = T),
                                          na_color = "bisque",
                                          reduction = "ft_umap", 
                                          features = "imputed_count_thres")
  p1 <- p1 + ggplot2::ggtitle(paste0(
    treatment, "\n", day_later, " growth potential of ", day_early, 
    " cells\n(UMAP of RNA fasttopics)\n(Log10-scale)")
  )
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_",
                                    treatment, "_imputation-ridge_ft-umap.png"),
                  p1, device = "png", width = 5, height = 5, units = "in")
  
  p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                          colors_use = list("red", "lightgray", "blue"),
                                          na_cutoff = quantile(seurat_object$imputed_count_thres, probs = 0.05, na.rm = T),
                                          na_color = "bisque",
                                          reduction = "umap", 
                                          features = "imputed_count_thres")
  p1 <- p1 + ggplot2::ggtitle(paste0(
    treatment, "\n", day_later, " growth potential of ", day_early, 
    " cells\n(UMAP of RNA fasttopics)\n(Log10-scale)")
  )
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_",
                                    treatment, "_imputation-ridge_umap.png"),
                  p1, device = "png", width = 5, height = 5, units = "in")
  
  ###########################
  
  print("Plotting lineage scatterplot")
  
  all(names(lineage_imputed_count) == names(lineage_future_count))
  
  lineage_imputed_count2 <- log10(lineage_imputed_count+1)
  lineage_future_count2 <- log10(lineage_future_count+1)
  
  labeling_vec <- rep(FALSE, length(lineage_imputed_count2))
  labeling_vec[intersect(which(lineage_imputed_count2 >= 1.5),
                         order(lineage_imputed_count2, decreasing = T)[1:10])] <- TRUE
  labeling_vec[intersect(which(lineage_future_count2 >= 1.5),
                         order(lineage_future_count2, decreasing = T)[1:10])] <- TRUE
  
  n <- length(lineage_imputed_count2)
  df <- data.frame(lineage_imputed_count = lineage_imputed_count2,
                   lineage_future_count = log10(lineage_future_count+1++ stats::runif(n, min = 0, max = 0.5)),
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
    " cells\n(Ridge for RNA fasttopics), (Log-scale)",
    "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2))
  ) +
    ggplot2::xlab("Observed lineage count (Log10, jittered)") + ggplot2::ylab("Predicted lineage count (Log10)")
  p1 <- p1 + Seurat::NoLegend() + ggplot2::coord_fixed()
  
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_", treatment, 
                                    "_imputation-ridge_lineage-level.png"),
                  p1, device = "png", width = 10, height = 10, units = "in")
  
  
  print("Saving")
  
  fit <- final_fit$fit
  
  save(cell_features,
       cell_imputed_score,
       fit, 
       lineage_imputed_count, 
       date_of_run, session_info,
       file = paste0("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_", treatment, "_from_day", day_early, "_postprocess.RData")
  )
  
  print("===========")
  
}

print("Done! :)")