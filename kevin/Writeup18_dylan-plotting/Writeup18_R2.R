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
  
  tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
  lineage_imputed_count <- all_data@misc[[paste0("fatepotential_", file)]]$lineage_imputed_count
  lineage_future_count <- tab_mat[names(lineage_imputed_count), day_later_full]
  
  cor_val <- stats::cor(lineage_imputed_count, lineage_future_count)
  print(paste0("Correlation: ", round(cor_val, 4)))
  print(paste0("Rsquared: ", round(cor_val^2, 4)))
  
  print("===")
}