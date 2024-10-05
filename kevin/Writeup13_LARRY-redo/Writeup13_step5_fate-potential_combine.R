rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup13/"
plot_folder <- "~/project/Multiome_fate/git/multiomeFate_analysis_kevin/fig/kevin/Writeup13/"

load(paste0(out_folder, "Writeup13_larry-dataset_step2_fasttopics.RData"))

file_vec <- c(Mono_d2_d4 = "Writeup13_Monocyte-4_from_day2_lineage-imputation.RData",
              Mono_d4_d6 = "Writeup13_Monocyte-6_from_day4_lineage-imputation.RData",
              Neu_d2_d4 = "Writeup13_Neutrophil-4_from_day2_lineage-imputation.RData",
              Neu_d4_d6 = "Writeup13_Neutrophil-6_from_day4_lineage-imputation.RData",
              Undiff_d2_d4 = "Writeup13_Undifferentiated-4_from_day2_lineage-imputation.RData",
              Undiff_d4_d6 = "Writeup13_Undifferentiated-6_from_day4_lineage-imputation.RData")

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
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup13_", filename, "_train-test-curve.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  # preparing the objects to be in seurat format
  fate_vec <- final_fit$cell_imputed_score
  cellnames <- Seurat::Cells(seurat_object)
  full_vec <- rep(NA, length(cellnames))
  names(full_vec) <- cellnames
  full_vec[names(fate_vec)] <- fate_vec
  seurat_object@meta.data[,paste0("fatepotential_", filename)] <- full_vec
  
  # put it into misc
  seurat_object@misc <- c(seurat_object@misc,
                     list(final_fit))
  names(seurat_object@misc)[length(seurat_object@misc)] <- paste0("fatepotential_", filename)
}

# save all the fate potential miscs and metadata
date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(seurat_object, date_of_run, session_info,
     file = paste0(out_folder, "Writeup13_combined_fatepotential.RData"))



