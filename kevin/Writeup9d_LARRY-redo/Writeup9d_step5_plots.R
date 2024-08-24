rm(list=ls())
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggtern)
library(multiomeFate)

data_path <- "~/project/Multiome_fate/out/kevin/Writeup9c/"
load(paste0(data_path, "Writeup9c_larry-dataset_step3_fasttopics.RData"))
date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

plot_path <- "~/project/Multiome_fate/out/figures/kevin/Writeup9d/"
plot_prefix <- paste0(plot_path, "Writeup9d_")

seurat_object_safe <- seurat_object

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
celltype_vec <- c("Monocyte", "Neutrophil", "Undifferentiated")
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
day_later <- "6"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  day_later <- treatment
  
  ###################
  
  data_path2 <- "~/project/Multiome_fate/out/kevin/Writeup9d/"
  load(paste0(data_path2, "Writeup9d_", treatment, "_from_day", day_early, "_lineage-imputation.RData"))
  class(fit_res) <- "lineage_cv"
  
  #####################
  
  print("Plotting train-test curve")
  plot1 <- multiomeFate:::plot_trainTest(cv_fit_list = fit_res,
                                         title_test = paste0(day_later, " growth potential of ", day_early, "\n(Testing)"),
                                         title_train = paste0(day_later, " growth potential of ", day_early, "\n(Training)"))
  ggplot2::ggsave(filename = paste0(plot_prefix, treatment, "_train-test-curve.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  ###########################
  
  print("Plotting ANOVA")
  
  plot1 <- multiomeFate:::plot_anova(seurat_object = seurat_object,
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
  
  plot1 <- multiomeFate:::plot_cellGrowthUmap(
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
  
  plot1 <- multiomeFate:::plot_lineageScatterplot(
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
       file = paste0(data_path2, "Writeup9d_", treatment, "_from_day", day_early, "_postprocess.RData")
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

df <- multiomeFate:::compute_entropy(cell_imputation_mat = cell_imputation_mat,
                                     later_timepoint = 6,
                                     seurat_object = seurat_object,
                                     variable_celltype = "Cell.type.annotation",
                                     variable_lineage = "assigned_lineage",
                                     variable_timepoint = "Time.point")

color_palette <- c("blue3", "coral2", "gray50")
names(color_palette) <- paste0(c("Monocyte", "Neutrophil", "Undifferentiated"))

for(color_by in c("celltype", "dominant_fate")){
  for(size_by in c("cellsize", "entropy")){
    
    color_by_name <- ifelse(color_by == "celltype", "celltype", "dominantFate")
    color_by_plot <- ifelse(color_by == "celltype", "current cell type", "future dominant fate")
    size_by_plot <- ifelse(size_by == "cellsize", "# predicted progenies", "entropy")
    
    aes_formula <- ggplot2::aes_string(x = "Monocyte", 
                                       y = "Neutrophil", 
                                       z = "Undifferentiated", 
                                       color = color_by, 
                                       size = size_by)
    # rearrange by color
    tmp <- lapply(rev(names(color_palette)), function(celltype){
      which(df[,color_by] == celltype)
    })
    df <- df[unlist(tmp),]
    
    # create empty plot
    pdf(file = paste0(plot_prefix, "d46_simplex_color-", color_by_name, "_size-", size_by, ".pdf"),
        width = 10, 
        height = 5, 
        onefile = TRUE)
    
    plot_list <- vector("list", 4)
    
    # plot of all
    # https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
    for(kk in 1:3){
      df_tmp <- df[which(df[,color_by] == names(color_palette)[kk]),]
      plot_list[[kk]] <- multiomeFate:::plot_simplex(
        df = df_tmp,
        aes_formula = aes_formula, 
        col_palette = color_palette,
        xlab = "Monocyte",
        ylab = "Neutrophil",
        zlab = "Undifferentiated",
        title = paste0("Day ", day_early, " cells predicting ", day_later, ", color by ", color_by_plot, 
                       "\nSize by ", size_by_plot,
                       "\nOnly ", names(color_palette)[kk])
      )
    }
    
    plot_list[[4]] <- multiomeFate:::plot_simplex(
      df = df,
      aes_formula = aes_formula, 
      col_palette = color_palette,
      xlab = "Monocyte",
      ylab = "Neutrophil",
      zlab = "Undifferentiated",
      title = paste0("Day ", day_early, " cells predicting ", day_later, ", color by ", color_by_plot, "\nSize by ", size_by_plot)
    )
    
    if(size_by == "entropy"){
      for(kk in 1:length(plot_list)){
        plot_list[[kk]] <- plot_list[[kk]] + ggplot2::scale_size_area(max_size = 1)
      }
    }
    
    for(kk in 1:length(plot_list)){
      print(plot_list[[kk]])
    }
    
    graphics.off()
  }
}


print("Done! :)")