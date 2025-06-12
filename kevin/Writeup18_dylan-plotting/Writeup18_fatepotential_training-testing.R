rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup18/"

all_data <- multiomeFate:::data_loader(which_files = c("lineage"))

file_vec <- c(CIS_d0_d10 = "Writeup10a_CIS-from-d0_fatepotential.RData",
              CIS_d10_w5 = "Writeup10a_CIS-from-d10_fatepotential.RData",
              COCL2_d0_d10 = "Writeup10a_COCL2-from-d0_fatepotential.RData",
              COCL2_d10_w5 = "Writeup10a_COCL2-from-d10_fatepotential.RData",
              DABTRAM_d0_d10 = "Writeup10a_DABTRAM-from-d0_fatepotential.RData",
              DABTRAM_d10_w5 = "Writeup10a_DABTRAM-from-d10_fatepotential.RData")

treatment_vec <- c("CIS", "CIS",
                   "COCL2", "COCL2",
                   "DABTRAM", "DABTRAM")
day_later_vec <- c("day10", "week5",
                   "day10", "week5",
                   "day10", "week5")

for(kk in 1:length(file_vec)){
  filepath <- file_vec[kk]
  filename <- names(file_vec)[kk]
  
  day_later <- day_later_vec[kk]
  treatment <- treatment_vec[kk]
  
  color_dataset <- all_data@misc$dataset_colors[[paste0(day_later, "_", treatment)]]
  
  print(paste0("Working on ", filename))
  # load the fate potential
  load(paste0(out_folder, filepath))
  
  # make the training-testing plots
  plot1 <- multiomeFate:::plot_trainTest(cv_fit_list = fit_res,
                                         fill_col = color_dataset,
                                         title_test = paste0(filename, " growth potential\n(Testing)"),
                                         title_train = paste0(filename, " growth potential\n(Training)"))
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup18_", filename, "_train-test-curve.png"),
                  plot1, width = 6, height = 3, units = "in")
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup18_", filename, "_train-test-curve.pdf"),
                  plot1, width = 5, height = 2.5, units = "in")
}
