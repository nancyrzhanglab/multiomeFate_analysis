rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup11/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup11/"

load(paste0(out_folder, "Writeup11_larry_fatepotential.RData"))

file_vec <- c("Monocyte_d2_d4", "Monocyte_d4_d6",
              "Neutrophil_d2_d4", "Neutrophil_d4_d6",
              "Undifferentiated_d2_d4", "Undifferentiated_d4_d6")
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
  
  cell_imputed_score <- seurat_obj@meta.data[,paste0("fatepotential_", file)]
  names(cell_imputed_score) <- Seurat::Cells(seurat_obj)
  
  ###########################
  
  print("Plotting violin plots")
  plot1 <- multiomeFate:::plot_anova(seurat_object = seurat_obj,
                                     cell_imputed_score = cell_imputed_score,
                                     assigned_lineage_variable = "assigned_lineage",
                                     time_celltype_variable = "time_celltype",
                                     day_later = day_later_full,
                                     ylab = paste0(day_later_vec, " growth potential"),
                                     ylim = c(max(stats::quantile(cell_imputed_score, na.rm = TRUE, prob = 0.05), -5), 
                                              NA))
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup11_", file, "_fatepotential-violinplot.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  
  plot1 <- multiomeFate:::plot_cellGrowthUmap(
    seurat_object = seurat_object,
    cell_imputed_score = cell_imputed_score,
    reduction = "python_X_emb",
    title = paste0(
      treatment, "\n", day_later, " growth potential of ", day_early, 
      " cells\n(UMAP of RNA fasttopics)\n(Log10-scale)")
  )
  
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup11_", file, "_fatepotential-umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  ###########################
  
  print("Plotting lineage scatterplot")
  
  lineage_imputed_count <- seurat_obj@misc[[paste0("fatepotential_", file)]]$lineage_imputed_count
  tab_mat <- table(seurat_obj$assigned_lineage, seurat_obj$time_celltype)
  lineage_future_count <- tab_mat[names(lineage_imputed_count), day_later_full]
  
  plot1 <- multiomeFate:::plot_lineageScatterplot(
    lineage_future_count = lineage_future_count,
    lineage_imputed_count = lineage_imputed_count,
    title = paste0(day_later, " growth potential of ", day_early, 
                   " cells\n(Ridge for RNA fasttopics), (Log-scale)")
  )
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup11_", file, "_fatepotential-lineage_prediction.png"),
                  plot1, device = "png", width = 10, height = 10, units = "in")
  
}

print("Done! :)")