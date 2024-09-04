rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"

all_data <- multiomeFate:::data_loader(which_files = c("fasttopics", "fatepotential"))

file_vec <- c("CIS_d0_d10", "CIS_d10_w5",
              "COCL2_d0_d10", "COCL2_d10_w5",
              "DABTRAM_d0_d10", "DABTRAM_d10_w5")
treatment_vec <- c("CIS", "CIS",
                   "COCL2", "COCL2",
                   "DABTRAM", "DABTRAM")
day_later_vec <- c("day10", "week5",
                   "day10", "week5",
                   "day10", "week5")
day_early_vec <- c("day0", "day10",
                   "day0", "day10",
                   "day0", "day10")

for(kk in 1:length(file_vec)){
  file <- file_vec[kk]
  day_later <- day_later_vec[kk]
  day_early <- day_early_vec[kk]
  treatment <- treatment_vec[kk]
  day_later_full <- paste0(day_later, "_", treatment)
  print(paste0("Working on ", day_later_full))
  
  ###########################
  
  print("Plotting violin plots")
  all_data2 <- all_data
  cell_imputed_score <- all_data2@meta.data[,paste0("fatepotential_", file)]
  keep_vec <- rep(FALSE, length(Seurat::Cells(all_data2)))
  keep_vec[!is.na(cell_imputed_score)] <- TRUE
  all_data2$keep <- keep_vec
  all_data2 <- subset(all_data2, keep == TRUE)
  
  cell_imputed_score <- all_data2@meta.data[,paste0("fatepotential_", file)]
  
  lineage_vec <- all_data2$assigned_lineage
  tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
  tab_vec <- table(lineage_vec)
  tab_vec <- tab_vec[tab_vec >= 2]
  later_size <- tab_mat[names(tab_vec), day_later_full]
  lineage_names <- names(sort(later_size, decreasing = T))[1:20]
  idx <- which(lineage_vec %in% lineage_names)
  
  df <- data.frame(lineage = lineage_vec[idx],
                   imputed_count = cell_imputed_score[idx])
  df_tmp <- df; df_tmp$lineage <- as.factor(df_tmp$lineage)
  anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df_tmp)
  df2 <- data.frame(lineage = "All",
                    imputed_count = cell_imputed_score)
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
  p1 <- p1 + ggplot2::ylab(paste0(day_later, " growth potential"))
  p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  p1 <- p1 + ggplot2::stat_summary(fun=max, geom="point", shape=10, size=5, color="blue")
  p1 <- p1 + ggplot2::ggtitle(paste0("ANOVA -Log10(pvalue)=", round(-log10(anova_res$p.value), 2), ", Lineage effect = ", lineage_effect, "%"))
  ggplot2::ggsave(filename = paste0(plot_folder, 
                                    "Writeup10a_",
                                    file, 
                                    "_fatepotential-violinplot.png"),
                  p1, device = "png", width = 6, height = 3, units = "in")
  
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  all_data2 <- all_data
  keep_vec <- rep(FALSE, length(Seurat::Cells(all_data2)))
  keep_vec[which(all_data2$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))] <- TRUE
  all_data2$keep <- keep_vec
  all_data2 <- subset(all_data2, keep == TRUE)
  
  cell_imputed_score <- all_data2@meta.data[,paste0("fatepotential_", file)]
  max_val <- stats::quantile(cell_imputed_score, probs = 0.99, na.rm = TRUE)
  all_data2$imputed_count_thres <- pmin(cell_imputed_score, max_val)
  
  p1 <- scCustomize::FeaturePlot_scCustom(all_data2, 
                                          colors_use = all_data@misc$fatepotential_colors,
                                          na_cutoff = quantile(all_data2$imputed_count_thres, 
                                                               probs = 0.05, 
                                                               na.rm = TRUE),
                                          na_color = all_data@misc$fatepotential_na_colors,
                                          reduction = paste0("ft.", treatment, ".umap"), 
                                          features = "imputed_count_thres")
  p1 <- p1 + ggplot2::ggtitle(paste0(
    treatment, "\n", day_later, " growth potential of ", day_early, 
    " cells\n(UMAP of RNA fasttopics)\n(Log10-scale)")
  )
  ggplot2::ggsave(filename = paste0(plot_folder, 
                                    "Writeup10a_",
                                    file, 
                                    "_fatepotential-umap.png"),
                  p1, device = "png", width = 5, height = 5, units = "in")
  
  ###########################
  
  print("Plotting lineage scatterplot")
  
  lineage_imputed_count <- all_data@misc[[paste0("fatepotential_", file)]]$lineage_imputed_count
  tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
  lineage_future_count <- tab_mat[names(lineage_imputed_count), day_later_full]
  
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
    " cells\n(Ridge for RNA fasttopics, ATAC PeakVI),",
    "\nCorrelation:", round(stats::cor(lineage_imputed_count2, lineage_future_count2), 2))
  ) +
    ggplot2::xlab("Observed lineage count (Log10, jittered)") + ggplot2::ylab("Predicted lineage count (Log10)")
  p1 <- p1 + Seurat::NoLegend() + ggplot2::coord_fixed()
  
  ggplot2::ggsave(filename = paste0(plot_folder, 
                                    "Writeup10a_",
                                    file, 
                                    "_fatepotential-lineage_prediction.png"),
                  p1, device = "png", width = 10, height = 10, units = "in")
  
}