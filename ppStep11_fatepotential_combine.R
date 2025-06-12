rm(list=ls())
library(Seurat)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
plot_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/figures/kevin/Writeup10a/"

all_data <- multiomeFate:::data_loader(which_files = c("lineage"))

file_vec <- c(CIS_d0_d10 = "Writeup10a_CIS-from-d0_fatepotential.RData",
              CIS_d10_w5 = "Writeup10a_CIS-from-d10_fatepotential.RData",
              COCL2_d0_d10 = "Writeup10a_COCL2-from-d0_fatepotential.RData",
              COCL2_d10_w5 = "Writeup10a_COCL2-from-d10_fatepotential.RData",
              DABTRAM_d0_d10 = "Writeup10a_DABTRAM-from-d0_fatepotential.RData",
              DABTRAM_d10_w5 = "Writeup10a_DABTRAM-from-d10_fatepotential.RData")

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
  ggplot2::ggsave(filename = paste0(plot_folder, "Writeup10a_", filename, "_train-test-curve.png"),
                  plot1, device = "png", width = 6, height = 3, units = "in")
  
  # preparing the objects to be in seurat format
  fate_vec <- final_fit$cell_imputed_score
  cellnames <- Seurat::Cells(all_data)
  full_vec <- rep(NA, length(cellnames))
  names(full_vec) <- cellnames
  full_vec[names(fate_vec)] <- fate_vec
  all_data@meta.data[,paste0("fatepotential_", filename)] <- full_vec
  
  # put it into misc
  all_data@misc <- c(all_data@misc,
                     list(final_fit))
  names(all_data@misc)[length(all_data@misc)] <- paste0("fatepotential_", filename)
}

# save all the fate potential miscs and metadata
date_of_run <- Sys.time()
session_info <- devtools::session_info()

all_data_fatepotential <- all_data@misc
all_data_fatepotential$dataset_colors <- NULL

save(all_data_fatepotential, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_fatepotential.RData"))

# Updating the Writeup10a_data_empty.RData
all_data_tmp <- all_data
load(paste0(out_folder, "Writeup10a_data_empty.RData"))
metadata_full <- all_data_tmp@meta.data
metadata_empty <- all_data@meta.data
cellnames_full <- Seurat::Cells(all_data_tmp)
cellnames_empty <- Seurat::Cells(all_data)

var_setdiff <- setdiff(colnames(metadata_full), colnames(metadata_empty))
for(variable in var_setdiff){
  vec <- rep(NA, length(cellnames_empty))
  names(vec) <- cellnames_empty
  vec[rownames(metadata_full)] <- metadata_full[,variable]
  metadata_empty[,variable] <- vec
}

all_data@meta.data <- metadata_empty

save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_data_empty.RData"))




