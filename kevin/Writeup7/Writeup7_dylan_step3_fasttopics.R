rm(list=ls())
library(Seurat)
library(fastTopics)

load("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step2_saver.RData")
set.seed(10)

treatment_vec <- sort(unique(all_data$OG_condition))
treatment_vec <- treatment_vec[treatment_vec != "naive"]

for(treatment in treatment_vec){
  print(paste0("Working on ", treatment))
  
  keep_vec <- rep(FALSE, ncol(all_data))
  keep_vec[which(all_data$dataset %in% c("naive", treatment))] <- TRUE
  all_data$keep <- keep_vec
  all_data2 <- subset(all_data, keep == TRUE)
  
  mat <- Seurat::GetAssayData(object = all_data2, assay = "RNA", slot = "counts")
  mat <- mat[Seurat::VariableFeatures(all_data2),]
  mat <- Matrix::t(mat)
  
  K <- 30
  set.seed(10)
  topic_res <- fastTopics::fit_topic_model(mat, k = K)
  
  #########
  
  topic_mat <- matrix(NA, nrow = ncol(all_data), ncol = ncol(topic_res$L))
  rownames(topic_mat) <- colnames(all_data)
  topic_res$L <- topic_res$L[rownames(topic_res$L) %in% colnames(all_data),]
  
  for(i in 1:nrow(topic_res$L)){
    topic_mat[rownames(topic_res$L)[i],] <- topic_res$L[i,]
  }
  colnames(topic_mat) <- paste0("fastTopic", treatment, "_", 1:ncol(topic_mat))
  colnames(topic_res$F) <- paste0("fastTopic", treatment, "_", 1:ncol(topic_mat))
  
  all_data[[paste0("fasttopic_", treatment)]] <- Seurat::CreateDimReducObject(embeddings = topic_mat, 
                                                                              loadings =  topic_res$F,
                                                                              assay = "RNA",
                                                                              key =  paste0("fastTopic", treatment, "_"))
}

print("Finished")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(all_data, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step3_fasttopics.RData")

print("Done! :)")


