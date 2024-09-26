rm(list=ls())
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggtern)
library(multiomeFate)

data_path <- "~/project/Multiome_fate/out/kevin/Writeup13/"
load(paste0(data_path, "Writeup13_larry-dataset_step2_fasttopics.RData"))
date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

plot_path <- "~/project/Multiome_fate/out/figures/kevin/Writeup13/"
plot_prefix <- paste0(plot_path, "Writeup13_")

cell_imputation_mat <- numeric(0)

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
celltype_vec <- c("Monocyte", "Neutrophil", "Undifferentiated")
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
day_later <- "6"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

for(treatment in treatment_vec){
  load(paste0(data_path, "Writeup13_", treatment, "_from_day", day_early, "_lineage-imputation.RData"))
  
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

