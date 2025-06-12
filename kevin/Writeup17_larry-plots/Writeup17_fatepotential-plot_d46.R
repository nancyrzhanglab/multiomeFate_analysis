rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup13/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17/"
code_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup17_larry-plots/"

load(paste0(out_folder, "Writeup13_combined_fatepotential.RData"))

file_vec <- c("Mono_d4_d6", "Neu_d4_d6", "Undiff_d4_d6")
treatment_vec <- c("Monocyte", "Neutrophil", "Undifferentiated")
day_later_vec <- rep("6", 3)
day_early_vec <- rep("4", 3)
color_terminal_vec <- c("#4472C4", "#D55E00", "#7F7F7F")

for(kk in 1:length(file_vec)){
  file <- file_vec[kk]
  day_later <- day_later_vec[kk]
  day_early <- day_early_vec[kk]
  treatment <- treatment_vec[kk]
  day_later_full <- paste0(treatment, "-", day_later)
  print(paste0("Working on ", day_later_full))
  
  ###########################
  
  print("Plotting UMAP with lineage imputation")
  
  seurat_object2 <- seurat_object
  keep_vec <- rep(FALSE, length(Seurat::Cells(seurat_object2)))
  keep_vec[seurat_object2$Time.point == as.numeric(day_early)] <- TRUE
  seurat_object2$keep <- keep_vec
  seurat_object2 <- subset(seurat_object2, keep == TRUE)
  
  embedding <- seurat_object2[["SPRING"]]@cell.embeddings
  
  cell_imputed_score <- seurat_object2@meta.data[,paste0("fatepotential_", file)]
  names(cell_imputed_score) <- Seurat::Cells(seurat_object2)

  min_val <- quantile(cell_imputed_score, probs = 0.1)
  max_val <- quantile(cell_imputed_score, probs = 0.9)
  
  color_palette <- grDevices::colorRampPalette(c("bisque", color_terminal_vec[kk]))(25)
  break_palette <- seq(min_val, max_val, length = 25)
  col_vec <- sapply(cell_imputed_score, function(x){
    color_palette[which.min(abs(x - break_palette))]
  })
  
  set.seed(10)
  order_idx <- sample(1:nrow(embedding))
  embedding <- embedding[order_idx,]
  cell_imputed_score <- cell_imputed_score[order_idx]
  col_vec <- col_vec[order_idx]
  
  plot1 <- multiomeFate:::plot_cellGrowthUmap(
    seurat_object = seurat_object2,
    cell_imputed_score = cell_imputed_score,
    colors_use = c("bisque", color_terminal_vec[kk]),
    na_color = "gray90",
    reduction = "SPRING",
    title = paste0(
      treatment, "\n", day_later, " growth potential of ", day_early, 
      " cells\n(UMAP of RNA fasttopics)\n(Log10-scale)")
  )
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17_", file,  "_fatepotential-umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  
  png(paste0(plot_folder, "Writeup17_", file,  "_fatepotential-umap_cleaned.png"),
      height = 1000, width = 1000, units = "px", res = 300)
  par(mar = rep(0.1,4))
  plot(embedding[,1],
       embedding[,2], 
       pch = 16,
       col = col_vec,
       xaxt = "n",
       yaxt = "n",
       bty = "n",
       xlab = "",
       ylab = "")
  graphics.off()
  
  ###########################
  
  print("Plotting lineage scatterplot")
  
  tab_mat <- table(seurat_object$assigned_lineage, seurat_object$time_celltype)
  lineage_imputed_count <- seurat_object@misc[[paste0("fatepotential_", file)]]$lineage_imputed_count
  lineage_future_count <- tab_mat[names(lineage_imputed_count), day_later_full]
  
  lineage_imputed_count2 <- log10(lineage_imputed_count+1)
  lineage_future_count2 <- log10(lineage_future_count+1)
  
  n <- length(lineage_imputed_count2)
  df <- data.frame(lineage_imputed_count = lineage_imputed_count2,
                   lineage_future_count = lineage_future_count2,
                   name = names(lineage_imputed_count2))
  
  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = lineage_future_count, 
                                            y = lineage_imputed_count))
  plot1 <- plot1 + ggplot2::geom_point()
  plot1 <- plot1 + Seurat::NoLegend() + ggplot2::coord_fixed()
  plot1 <- plot1 + ggplot2::labs(
    x = "Obs. lineage size (Log10)",
    y = "Pred. lineage size (Log10)"
  ) + Seurat::NoLegend() 
  plot1
  
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17_",  file, "_fatepotential-lineage_prediction_cleaned.png"),
                  plot1, device = "png", width = 2.5, height = 3, units = "in")
  
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Corr: ", round(stats::cor(lineage_imputed_count2, lineage_future_count2) , 2)))
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup17_",  file, "_fatepotential-lineage_prediction.png"),
                  plot1, device = "png", width = 4, height = 4, units = "in")
  
}