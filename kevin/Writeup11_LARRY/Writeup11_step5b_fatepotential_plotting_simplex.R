rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup11/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup11/"

load(paste0(out_folder, "Writeup11_larry_fatepotential.RData"))

# construct the matrix of fate potentials

time_early_vec <- c(2, 4)
time_later_vec <- c(4, 6)

for(kk in 1:length(time_early_vec)){
  time_early <- time_early_vec[kk]
  time_later <- time_later_vec[kk]
  
  cell_idx <- which(seurat_obj$time_info == time_early)
  
  colnames_vec <- colnames(seurat_obj@meta.data)
  colnames_fate_potential <- colnames_vec[grep(paste0("fatepotential_.*d", time_early, "_.*"), colnames_vec)]
  cell_imputation_mat <- seurat_obj@meta.data[cell_idx,colnames_fate_potential]
  
  stopifnot(all(!is.na(cell_imputation_mat)))
  
  celltype_vec <- sapply(colnames(cell_imputation_mat), function(x){
    strsplit(x, split = "_")[[1]][2]
  })
  colnames(cell_imputation_mat) <- celltype_vec
  
  df <- multiomeFate:::compute_entropy(cell_imputation_mat = cell_imputation_mat,
                                       later_timepoint = time_later,
                                       seurat_object = seurat_obj,
                                       variable_celltype = "state_info",
                                       variable_lineage = "assigned_lineage",
                                       variable_timepoint = "time_info")
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
          title = paste0("Day ", time_early, " cells predicting ", time_later, ", color by ", color_by_plot, 
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
        title = paste0("Day ", time_early, " cells predicting ", time_later, ", color by ", color_by_plot, "\nSize by ", size_by_plot)
      )
      
      if(size_by == "entropy"){
        for(kk in 1:length(plot_list)){
          plot_list[[kk]] <- plot_list[[kk]] + ggplot2::scale_size_area(max_size = 1)
        }
      }
      
      # create empty plot
      pdf(file = paste0(plot_folder, "Writeup11_time-", time_early, "-to-", time_later, "_simplex_color-", color_by_name, "_size-", size_by, ".pdf"),
          width = 10, 
          height = 5, 
          onefile = TRUE)
      
      for(kk in 1:length(plot_list)){
        print(plot_list[[kk]])
      }
      
      graphics.off()
    }
  }
}

print("Done! :)")