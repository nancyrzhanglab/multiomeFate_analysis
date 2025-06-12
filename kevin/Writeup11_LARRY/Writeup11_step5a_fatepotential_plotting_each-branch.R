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

#################

var_vec <- c(time_info = TRUE,
             state_info = TRUE,
             time_celltype = TRUE)

pdf(paste0(plot_folder, "Writeup11_larry_umap-covariates.pdf"),
    onefile = T, width = 6, height = 5)

for(kk in 1:length(var_vec)){
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = "python_X_emb",
                           group.by = names(var_vec)[kk],
                           raster = TRUE,
                           label = var_vec[kk])
  print(plot1)
}

dev.off()

#################

# Gini plots

unique_lineages <- sort(unique(seurat_obj$assigned_lineage))

treatment_vec <- unique(seurat_obj$time_celltype)
lineage_idx_list <- lapply(unique_lineages, function(lineage){
  which(seurat_obj$assigned_lineage == lineage)
})

pdf(paste0(plot_folder, "Writeup11_gini-index_by-time-celltype.pdf"),
    onefile = T, width = 5, height = 5)

for(treatment in treatment_vec){
  # count
  cell_treatment_idx <- which(seurat_obj$time_celltype == treatment)
  
  vec <- sapply(lineage_idx_list, function(lineage_idx){
    length(intersect(cell_treatment_idx,
                     lineage_idx))
  })
  
  gini_val <- dineq::gini.wtd(vec)
  
  plot(x = seq(0, 1, length.out = length(vec)), 
       y = cumsum(sort(vec, decreasing = FALSE))/sum(vec),
       pch = 16,
       asp = TRUE,
       xlab = "Cumulative share of lineages (smallest to largest)",
       ylab = "Cumulative share of cells",
       main = paste0("Gini index for ", treatment, ": ", round(gini_val, 2)))
  lines(c(0,1), 
        c(0,1), 
        col = "red",
        lwd = 2, 
        lty = 2)
}

dev.off()

####

treatment_vec <- sort(unique(seurat_obj$time_info))
pdf(paste0(plot_folder, "Writeup11_gini-index_by-time.pdf"),
    onefile = T, width = 5, height = 5)

for(treatment in treatment_vec){
  # count
  cell_treatment_idx <- which(seurat_obj$time_info == treatment)
  
  vec <- sapply(lineage_idx_list, function(lineage_idx){
    length(intersect(cell_treatment_idx,
                     lineage_idx))
  })
  
  gini_val <- dineq::gini.wtd(vec)
  
  plot(x = seq(0, 1, length.out = length(vec)), 
       y = cumsum(sort(vec, decreasing = FALSE))/sum(vec),
       pch = 16,
       asp = TRUE,
       xlab = "Cumulative share of lineages (smallest to largest)",
       ylab = "Cumulative share of cells",
       main = paste0("Gini index for ", treatment, ": ", round(gini_val, 2)))
  lines(c(0,1), 
        c(0,1), 
        col = "red",
        lwd = 2, 
        lty = 2)
}

dev.off()

print("Done! :)")