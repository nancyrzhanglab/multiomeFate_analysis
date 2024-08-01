rm(list=ls())
library(Seurat)
library(multiomeFate)
data_path <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup9c/"
load(paste0(data_path, "Writeup9c_larry-dataset_step3_fasttopics.RData"))
date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

plot_path <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/figures/Writeup9d/"
plot_prefix <- paste0(plot_path, "Writeup9d_")

seurat_object_safe <- seurat_object

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  day_later <- treatment
  
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
  
  plot1 <- Seurat::DimPlot(seurat_object, reduction = "ft_umap",
                           group.by = "time_celltype", pt.size = .3, label = T)
  ggplot2::ggsave(filename = paste0(plot_prefix, treatment, "_umap.png"),
                  plot1, device = "png", width = 6, height = 4, units = "in")
  
  ###################
  
  load(paste0(data_path, "Writeup9c_", treatment, "_from_day", day_early, "_lineage-imputation.RData"))
  class(fit_res) <- "lineage_cv"
  
  #####################
  
  print("Plotting train-test curve")
  plot1 <- plot_trainTest(cv_fit_list = fit_res,
                          title_test = paste0(day_later, " growth potential of ", day_early, "\n(Testing)"),
                          title_train = paste0(day_later, " growth potential of ", day_early, "\n(Training)"))
  ggplot2::ggsave(filename = paste0(plot_prefix, treatment, "_train-test-curve.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  ###########################
  plot1 <- plot_anova(seurat_object = seurat_object,
                      cell_imputed_score = cell_imputed_score,
                      assigned_lineage_variable = "assigned_lineage",
                      time_celltype_variable = "time_celltype",
                      day_later = day_later,
                      ylab = "Week5 growth potential")
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9c/Writeup9c_",
                                    treatment, "_lineage-growthpotential-violinplot.png"),
                  p1, device = "png", width = 6, height = 3, units = "in")
  
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  
  cell_imputed_score_full <- rep(NA, ncol(seurat_object))
  names(cell_imputed_score_full) <- colnames(seurat_object)
  cell_imputed_score_full[names(cell_imputed_score)] <- cell_imputed_score
  
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
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9c/Writeup9c_",
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
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9c/Writeup9c_",
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
  
  ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9c/Writeup9c_", treatment, 
                                    "_imputation-ridge_lineage-level.png"),
                  p1, device = "png", width = 10, height = 10, units = "in")
  
  
  print("Saving")
  
  fit <- final_fit
  
  save(cell_features,
       cell_imputed_score,
       fit, 
       lineage_imputed_count, 
       date_of_run, session_info,
       file = paste0("~/project/Multiome_fate/out/kevin/Writeup9c/Writeup9c_", treatment, "_from_day", day_early, "_postprocess.RData")
  )
  
  print("===========")
  
}

print("Done! :)")