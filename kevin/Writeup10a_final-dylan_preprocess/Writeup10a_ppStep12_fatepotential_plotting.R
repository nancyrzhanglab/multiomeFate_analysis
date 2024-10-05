rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup10a/"

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
  
  cell_imputed_score <- all_data@meta.data[,paste0("fatepotential_", file)]
  names(cell_imputed_score) <- Seurat::Cells(all_data)
  
  plot1 <- multiomeFate:::plot_anova(seurat_object = all_data,
                                     cell_imputed_score = cell_imputed_score,
                                     assigned_lineage_variable = "assigned_lineage",
                                     time_celltype_variable = "dataset",
                                     day_later = day_later_full,
                                     ylab = paste0(day_later_vec, " growth potential"),
                                     ylim = c(max(stats::quantile(cell_imputed_score, na.rm = TRUE, prob = 0.05), -5), 
                                              NA))
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_", file, "_fatepotential-violinplot.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  
  all_data2 <- all_data
  keep_vec <- rep(FALSE, length(Seurat::Cells(all_data2)))
  keep_vec[which(all_data2$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))] <- TRUE
  all_data2$keep <- keep_vec
  all_data2 <- subset(all_data2, keep == TRUE)
  
  plot1 <- multiomeFate:::plot_cellGrowthUmap(
    seurat_object = all_data2,
    cell_imputed_score = cell_imputed_score,
    colors_use = all_data2@misc$fatepotential_colors,
    na_color = all_data2@misc$fatepotential_na_colors,
    reduction = paste0("ft.", treatment, ".umap"),
    title = paste0(
      treatment, "\n", day_later, " growth potential of ", day_early, 
      " cells\n(UMAP of RNA fasttopics)\n(Log10-scale)")
  )
  
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_", file,  "_fatepotential-umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  ###########################
  
  print("Plotting lineage scatterplot")
  
  
  tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
  lineage_imputed_count <- all_data@misc[[paste0("fatepotential_", file)]]$lineage_imputed_count
  lineage_future_count <- tab_mat[names(lineage_imputed_count), day_later_full]

  plot1 <- multiomeFate:::plot_lineageScatterplot(
    lineage_future_count = lineage_future_count,
    lineage_imputed_count = lineage_imputed_count,
    title = paste0(
      treatment, " ", day_later, " growth potential of ", day_early, 
      " cells\n(Ridge for RNA fasttopics, ATAC PeakVI), (Log-scale)")
  )
  
 
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_",  file, "_fatepotential-lineage_prediction.png"),
                  plot1, device = "png", width = 10, height = 10, units = "in")
  
}