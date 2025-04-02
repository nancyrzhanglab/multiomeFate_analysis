rm(list=ls())
library(Seurat)
library(ggplot2)
library(reshape2)

out_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/out/Writeup13/"
plot_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/fig/kevin/Writeup17/"
code_folder <- "/Users/kevinlin/Library/CloudStorage/Dropbox/Collaboration-and-People/Nancy/multiomeFate/git/multiomeFate_analysis/kevin/Writeup17_larry-plots/"

load(paste0(out_folder, "Writeup13_larry-dataset.RData"))

unique_lineages <- sort(unique(seurat_object$assigned_lineage))

lineage_idx_list <- lapply(unique_lineages, function(lineage){
  which(seurat_object$assigned_lineage == lineage)
})

######

# count
cell_treatment_idx <- which(seurat_object$Time.point == 6)

vec <- sapply(lineage_idx_list, function(lineage_idx){
  length(intersect(cell_treatment_idx,
                   lineage_idx))
})

gini <- dineq::gini.wtd(vec)

png(paste0(plot_folder, "Writeup17_gini.png"),
    height = 1200, width = 1200, units = "px", res = 300)
plot(NA,
     xlim = c(0,1),
     ylim = c(0,1),
     asp = TRUE,
     xlab = "Cumulative share of lineages (smallest to largest)",
     ylab = "Cumulative share of cells",
     main = paste0("Gini index (Day 6): ",
                   paste0(round(gini, 2), collapse = ", ")))
lines(x = seq(0, 1, length.out = length(vec)), 
      y = cumsum(sort(vec, decreasing = FALSE))/sum(vec),
      lwd = 3)
graphics.off()

png(paste0(plot_folder, "Writeup18_gini_cleaned.png"),
    height = 350, width = 350, units = "px", res = 300)
par(mar = c(0.5,0.5,0.1,0.1))
plot(NA,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "",
     ylab = "",
     asp = TRUE,
     xaxt = "n",
     yaxt = "n",
     bty = "n")
axis(1); axis(2)
lines(c(0,1), c(0,1), lty = 2, col = "coral", lwd = 3)

lines(x = seq(0, 1, length.out = length(vec)), 
      y = cumsum(sort(vec, decreasing = FALSE))/sum(vec),
      lwd = 8)
graphics.off()