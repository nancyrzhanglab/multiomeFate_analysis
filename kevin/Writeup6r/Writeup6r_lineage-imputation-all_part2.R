rm(list=ls())
library(Seurat)
library(Signac)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
day_early_vec <- c("day0", "day10")
day_later_vec <- c("day10", "week5")

#######################

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
 
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
  ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6r/Writeup6r_", treatment, "_umap.png"),
                  plot1, device = "png", width = 5, height = 4, units = "in")
  
  for(ii in 1:2){
    day_early <- day_early_vec[ii]
    day_later <- day_later_vec[ii]
    
    if(day_early == "day0") {
      day_early_full <- day_early
    } else {
      day_early_full <- paste0(day_early, "_", treatment)
    }
    day_later_full <- paste0(day_later, "_", treatment)
    
    print(paste0("Working on ", day_early))
    
    load(paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation.RData"))
    all_data2 <- all_data
    
    ###################
    
    train_mat <- sapply(cv_fit_list, function(x){
      x$train_loglik
    })
    
    test_mat <- sapply(cv_fit_list, function(x){
      x$test_loglik
    })
    
    train_mean <- Matrix::rowMeans(train_mat)
    test_mean <- Matrix::rowMeans(test_mat)
    train_sd <- apply(train_mat, 1, stats::sd)
    test_sd <- apply(test_mat, 1, stats::sd)
    
    lambda_sequence <- cv_fit_list[[1]]$train_fit$lambda_sequence
    lambda <- lambda_sequence[which.min(test_mean)]
    lambda
    
    final_fit <- multiomeFate:::lineage_imputation(
      cell_features = cell_features,
      cell_lineage = cell_lineage,
      coefficient_initial_list = cv_fit_list[[1]]$train_fit$fit_list[[which.min(test_mean)]]$coefficient_vec,
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
    
    png(paste0("../../../../out/figures/kevin/Writeup6r/Writeup6r_",
               treatment, "-", day_early, "_train-test-curve.png"),
        height = 1500, width = 3000, units = "px", res = 300)
    par(mfrow = c(1,2), mar = c(4,4,4,4))
    plot(lambda_sequence2, train_mean, 
         log = "x", type = "n",
         main = paste0(treatment, " ", day_later, " growth potential of ", day_early),
         xaxt = "n",
         xlab = "Lambda+1 (Log-scale tickmarks)", 
         ylab = "Negative loglikelihood (Training)",
         ylim = range(c(train_mean + train_sd, rev(train_mean - train_sd))))
    axis(side = 1, 
         at = xaxt_values, 
         labels = FALSE)
    axis(side = 1, 
         lwd = 2,
         at = 10^(0:max_base10), 
         labels = sapply(10^(0:max_base10), function(x){
           as.character(format(x, scientific = FALSE))
         }))
    polygon(x = c(lambda_sequence2, rev(lambda_sequence2)),
            y = c(train_mean + train_sd, rev(train_mean - train_sd)),
            col = "gray", 
            border = "black")
    points(lambda_sequence2, train_mean, pch = 16)
    lines(lambda_sequence2, train_mean, lwd = 2)
    
    plot(lambda_sequence2, test_mean, 
         log = "x", type = "n",
         main = paste0(treatment, " ", day_later, " growth potential of ", day_early,
                       "\nLambda = ", round(lambda,3)),
         xaxt = "n",
         xlab = "Lambda+1 (Log-scale tickmarks)", 
         ylab = "Negative loglikelihood (Testing)",
         ylim = range(c(test_mean + test_sd, rev(test_mean - test_sd))))
    axis(side = 1, 
         at = xaxt_values, 
         labels = FALSE)
    axis(side = 1, 
         lwd = 2,
         at = 10^(0:max_base10), 
         labels = sapply(10^(0:max_base10), function(x){
           as.character(format(x, scientific = FALSE))
         }))
    polygon(x = c(lambda_sequence2, rev(lambda_sequence2)),
            y = c(test_mean + test_sd, rev(test_mean - test_sd)),
            col = "gray", 
            border = "black")
    points(lambda_sequence2, test_mean, pch = 16)
    lines(lambda_sequence2, test_mean, lwd = 2)
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
    
    lineage_future_count <- tab_mat[uniq_lineage,day_later_full]
    
    ###########################
    
    print("Plotting UMAP with lineage imputation")
    
    cell_imputed_score_full <- rep(NA, ncol(all_data2))
    names(cell_imputed_score_full) <- colnames(all_data2)
    cell_imputed_score_full[rownames(cell_features)] <- cell_imputed_score
    
    all_data2$imputed_count <- cell_imputed_score_full
    
    p1 <- scCustomize::FeaturePlot_scCustom(all_data2, 
                                            colors_use = list("red", "lightgray", "blue"),
                                            na_cutoff = quantile(all_data2$imputed_count, probs = 0.05, na.rm = T),
                                            na_color = "bisque",
                                            reduction = "umap", 
                                            features = "imputed_count")
    p1 <- p1 + ggplot2::ggtitle(paste0(
      treatment, "\n", day_later, " growth potential of ", day_early, 
      " cells\n(RNA fasttopics, ATAC PeakVI)\n(Log-scale)")
    )
    ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6r/Writeup6r_",
                                      treatment, "-", day_early, "_imputation-ridge_umap.png"),
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
    p1 <- p1 + Seurat::NoLegend() + ggplot2::coord_fixed()
    
    ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6r/Writeup6r_", treatment, "-", day_early, "_imputation-ridge_lineage-level.png"),
                    p1, device = "png", width = 10, height = 10, units = "in")
    
    
    print("Saving")
    
    fit <- final_fit$fit
    
    save(cell_imputed_score,
         fit, 
         lineage_imputed_count, 
         date_of_run, session_info,
         file = paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation_postprocess.RData"))
    
    print("===========")
  }
}

print("Done! :)")

