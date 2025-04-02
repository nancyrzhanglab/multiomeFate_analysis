rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"

#######

all_data <- multiomeFate:::data_loader(which_files = c("fasttopics"))

unique_lineages <- sort(unique(all_data$assigned_lineage))

lineage_idx_list <- lapply(unique_lineages, function(lineage){
  which(all_data$assigned_lineage == lineage)
})

dataset_colors <- all_data@misc$dataset_colors

treatment_list <- list(CIS = c("day0", "day10_CIS", "week5_CIS"),
                       COCL2 = c("day0", "day10_COCL2", "week5_COCL2"),
                       DABTRAM = c("day0", "day10_DABTRAM", "week5_DABTRAM"))

for(kk in 1:length(treatment_list)){
  print(paste0("Working on: ", names(treatment_list)[kk]))
  treatment_vec <- treatment_list[[kk]]
  
  # count
  cell_treatment_idx_list <- lapply(treatment_vec, function(treatment){
    which(all_data$dataset == treatment)
  })
  
  vec_list <- lapply(cell_treatment_idx_list, function(cell_treatment_idx){
    sapply(lineage_idx_list, function(lineage_idx){
      length(intersect(cell_treatment_idx,
                       lineage_idx))
    })
  })
  
  gini_vec <- sapply(vec_list, function(vec){
    dineq::gini.wtd(vec)
  })
  
  png(paste0(plot_folder, "Writeup18_gini_", names(treatment_list)[kk], ".png"),
      height = 1200, width = 1200, units = "px", res = 300)
  plot(NA,
       xlim = c(0,1),
       ylim = c(0,1),
       asp = TRUE,
       xlab = "Cumulative share of lineages (smallest to largest)",
       ylab = "Cumulative share of cells",
       main = paste0("Gini index for ", names(treatment_list)[kk], ": ",
                     paste0(round(gini_vec, 2), collapse = ", ")))
  for(i in 1:length(treatment_vec)){
    lines(x = seq(0, 1, length.out = length(vec_list[[i]])), 
          y = cumsum(sort(vec_list[[i]], decreasing = FALSE))/sum(vec_list[[i]]),
          col = dataset_colors[treatment_vec[i]],
          lwd = 3)
  }
  graphics.off()
  
  png(paste0(plot_folder, "Writeup18_gini_", names(treatment_list)[kk], "_cleaned.png"),
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
  
  for(i in 1:length(treatment_vec)){
    lines(x = seq(0, 1, length.out = length(vec_list[[i]])), 
          y = cumsum(sort(vec_list[[i]], decreasing = FALSE))/sum(vec_list[[i]]),
          col = dataset_colors[treatment_vec[i]],
          lwd = ifelse(i == 3, 8, 4))
  }
  graphics.off()
}