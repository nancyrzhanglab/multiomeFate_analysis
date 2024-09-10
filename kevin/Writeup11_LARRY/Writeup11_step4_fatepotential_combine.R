rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup11/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup11/"
load(paste0(out_folder, "Writeup11_larry_seurat_fasttopics.RData"))

file_vec <- c(Monocyte_d2_d4 = "Writeup11_fatepotential_for_Monocyte-4.RData",
              Monocyte_d4_d6 = "Writeup11_fatepotential_for_Monocyte-6.RData",
              Neutrophil_d2_d4 = "Writeup11_fatepotential_for_Neutrophil-4.RData",
              Neutrophil_d4_d6 = "Writeup11_fatepotential_for_Neutrophil-6.RData",
              Undifferentiated_d2_d4 = "Writeup11_fatepotential_for_Undifferentiated-4.RData",
              Undifferentiated_d4_d6 = "Writeup11_fatepotential_for_Undifferentiated-6.RData")

for(kk in 1:length(file_vec)){
  filepath <- file_vec[kk]
  filename <- names(file_vec)[kk]
  
  print(paste0("Working on ", filename))
  # load the fate potential
  load(paste0(out_folder, filepath))
  
  # make the training-testing plots
  plot1 <- multiomeFate:::plot_trainTest(cv_fit_list = fit_res,
                                         title_test = paste0(filename, " growth potential\n(Testing)"),
                                         title_train = paste0(filename, " growth potential\n(Training)"))
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup11_", filename, "_train-test-curve.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  # preparing the objects to be in seurat format
  fate_vec <- final_fit$cell_imputed_score
  cellnames <- Seurat::Cells(seurat_obj)
  full_vec <- rep(NA, length(cellnames))
  names(full_vec) <- cellnames
  full_vec[names(fate_vec)] <- fate_vec
  seurat_obj@meta.data[,paste0("fatepotential_", filename)] <- full_vec
  
  # put it into misc
  seurat_obj@misc <- c(seurat_obj@misc,
                     list(final_fit))
  names(seurat_obj@misc)[length(seurat_obj@misc)] <- paste0("fatepotential_", filename)
}

# save all the fate potential miscs and metadata
date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(seurat_obj, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11_larry_fatepotential.RData"))
