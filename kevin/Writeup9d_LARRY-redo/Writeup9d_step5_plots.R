rm(list=ls())
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggtern)
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
celltype_vec <- c("Monocyte", "Neutrophil", "Undifferentiated")
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
  
  print("Plotting ANOVA")
  
  plot1 <- plot_anova(seurat_object = seurat_object,
                      cell_imputed_score = cell_imputed_score,
                      assigned_lineage_variable = "assigned_lineage",
                      time_celltype_variable = "time_celltype",
                      day_later = day_later,
                      ylab = "Week5 growth potential",
                      ylim = c(-5, NA))
  ggplot2::ggsave(filename = paste0(plot_prefix, treatment, "_lineage-growthpotential-violinplot.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  
  plot1 <- plot_cellGrowthUmap(
    seurat_object = seurat_object,
    cell_imputed_score = cell_imputed_score,
    reduction = "umap",
    title = paste0(
      treatment, "\n", day_later, " growth potential of ", day_early, 
      " cells\n(UMAP of RNA)\n(Log10-scale)")
  )
  ggplot2::ggsave(filename = paste0(plot_prefix, treatment, "_imputation-ridge_umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  ###########################
  
  print("Plotting lineage scatterplot")
  
  plot1 <- plot_lineageScatterplot(
    lineage_future_count = lineage_future_count,
    lineage_imputed_count = lineage_imputed_count,
    title = paste0(day_later, " growth potential of ", day_early, 
                   " cells\n(Ridge for RNA fasttopics), (Log-scale)")
  )
  ggplot2::ggsave(filename = paste0(plot_prefix, treatment, 
                                    "_imputation-ridge_lineage-level.png"),
                  plot1, device = "png", width = 10, height = 10, units = "in")
  
  ############################
  
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

########

cell_imputation_mat <- numeric(0)

for(treatment in treatment_vec){
  load(paste0(data_path, "Writeup9c_", treatment, "_from_day", day_early, "_postprocess.RData"))
  
  cell_imputation_mat <- cbind(cell_imputation_mat, cell_imputed_score)
}
colnames(cell_imputation_mat) <- celltype_vec

df <- compute_entropy(cell_imputation_mat = cell_imputation_mat,
                      later_timepoint = 6,
                      seurat_object = seurat_object,
                      variable_celltype = "Cell.type.annotation",
                      variable_lineage = "assigned_lineage",
                      variable_timepoint = "Time.point")

color_palette <- c("blue3", "coral2", "gray50")
names(color_palette) <- paste0(c("Monocyte", "Neutrophil", "Undifferentiated"))
df <- df[order(df$entropy, decreasing = TRUE),]

plot1 <- plot_simplex(
  df = df,
  aes_formula = ggplot2::aes(x = Monocyte, 
                             y = Neutrophil, 
                             z = Undifferentiated, 
                             color = celltype, 
                             size = cellsize), 
  col_palette = color_palette,
  xlab = "Monocyte",
  ylab = "Neutrophil",
  zlab = "Undifferentiated",
  title = paste0("Entropy of observed lineage for day ", day_early)
)

ggtern::ggsave(filename = paste0(plot_prefix, "percent-fate_", day_early, "_cell-pairs_by-entropy.png"),
               plot1, 
               device = "png", 
               width = 10, 
               height = 5, 
               units = "in")

print("Done! :)")