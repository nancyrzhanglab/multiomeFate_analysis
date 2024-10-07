rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup13/"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup13/"

load(paste0(out_folder, "Writeup13_combined_fatepotential.RData"))

file_vec <- c("Mono_d2_d4", "Mono_d4_d6",
              "Neu_d2_d4", "Neu_d4_d6",
              "Undiff_d2_d4", "Undiff_d4_d6")
treatment_vec <- c("Monocyte", "Monocyte",
                   "Neutrophil", "Neutrophil",
                   "Undifferentiated", "Undifferentiated")
day_later_vec <- c("4", "6",
                   "4", "6",
                   "4", "6")
day_early_vec <- c("2", "4",
                   "2", "4",
                   "2", "4")

for(kk in 1:length(file_vec)){
  file <- file_vec[kk]
  day_later <- day_later_vec[kk]
  day_early <- day_early_vec[kk]
  treatment <- treatment_vec[kk]
  day_later_full <- paste0(treatment, "-", day_later)
  print(paste0("Working on ", day_later_full))
  
  ###########################
  
  print("Plotting violin plots")
  
  cell_imputed_score <- seurat_object@meta.data[,paste0("fatepotential_", file)]
  names(cell_imputed_score) <- Seurat::Cells(seurat_object)
  
  plot1 <- multiomeFate:::plot_anova(seurat_object = seurat_object,
                                     cell_imputed_score = cell_imputed_score,
                                     assigned_lineage_variable = "assigned_lineage",
                                     time_celltype_variable = "time_celltype",
                                     day_later = day_later_full,
                                     ylab = paste0(day_later_vec, " growth potential"),
                                     ylim = c(max(stats::quantile(cell_imputed_score, na.rm = TRUE, prob = 0.05), -5), 
                                              NA))
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup13_", file, "_fatepotential-violinplot.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  
  seurat_object2 <- seurat_object
  keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_object2)))
  keep_vec[seurat_object2$Time.point == as.numeric(day_early)] <- TRUE
  seurat_object2$keep <- keep_vec
  seurat_object2 <- subset(seurat_object2, keep == TRUE)
  
  plot1 <- multiomeFate:::plot_cellGrowthUmap(
    seurat_object = seurat_object2,
    cell_imputed_score = cell_imputed_score,
    colors_use = c("red", "bisque", "blue"),
    na_color = "gray90",
    reduction = "SPRING",
    title = paste0(
      treatment, "\n", day_later, " growth potential of ", day_early, 
      " cells\n(UMAP of RNA fasttopics)\n(Log10-scale)")
  )
  
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup13_", file,  "_fatepotential-umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  ###########################
  
  print("Plotting lineage scatterplot")
  
  tab_mat <- table(seurat_object$assigned_lineage, seurat_object$time_celltype)
  lineage_imputed_count <- seurat_object@misc[[paste0("fatepotential_", file)]]$lineage_imputed_count
  lineage_future_count <- tab_mat[names(lineage_imputed_count), day_later_full]
  
  plot1 <- multiomeFate:::plot_lineageScatterplot(
    lineage_future_count = lineage_future_count,
    lineage_imputed_count = lineage_imputed_count,
    title = paste0(
      treatment, " ", day_later, " growth potential of ", day_early, 
      " cells\n(Ridge for RNA fasttopics), (Log-scale)")
  )
  
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup13_",  file, "_fatepotential-lineage_prediction.png"),
                  plot1, device = "png", width = 10, height = 10, units = "in")
  
}