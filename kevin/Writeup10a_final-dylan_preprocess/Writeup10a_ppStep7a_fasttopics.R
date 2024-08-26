rm(list=ls())
library(Seurat)
library(Signac)
library(fastTopics)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment_vec <- c("CIS", "COCL2", "DABTRAM")

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  keep_vec <- rep(FALSE, ncol(all_data))
  keep_vec[which(all_data$dataset %in% c("day0", 
                                         paste0("day10_", treatment), 
                                         paste0("week5_", treatment)))] <- TRUE
  all_data$keep <- keep_vec
  all_data2 <- subset(all_data, keep == TRUE)
  
  mat <- SeuratObject::LayerData(all_data2,
                                 assay = "RNA",
                                 features = Seurat::VariableFeatures(all_data2, assay = "RNA"),
                                 layer = "counts")
  mat <- Matrix::t(mat)
  
  print("Computing fastTopics")
  K <- 30
  set.seed(10)
  topic_res <- fastTopics::fit_topic_model(mat, k = K)
  
  print("Saving")
  save(topic_res, date_of_run, session_info,
       file = paste0(out_folder, "Writeup10a_ppStep7a_fasttopic_", treatment, ".RData"))
}

print("Done! :)")
