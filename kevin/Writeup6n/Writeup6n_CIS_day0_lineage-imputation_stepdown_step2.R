rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

treatment <- "CIS"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_day0_lineage-imputation_stepdown.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

loocv_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("day10_", treatment)] >= 20),
                                              which(tab_mat[,"day0"] >= 1))]
loocv_len <- length(loocv_lineages)
loocv_idx_list <- lapply(loocv_lineages, function(lineage){
  which(cell_lineage == lineage)
})

#########################

len <- length(coefficient_list_list)
loocv_mat <- matrix(NA, nrow = loocv_len, ncol = len)

for(iteration in 1:len){
  print(paste0("On iteration ", iteration))
  var_names <- coefficient_list_list[[iteration]]$var_next_iteration
  
  cell_features <- cell_features_full[,var_names,drop=F]
  p <- ncol(cell_features)
  tmp <- quantile(abs(cell_features[,-which(colnames(cell_features) == "Intercept"),drop=F]), probs = 0.95)
  coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
  coefficient_initial <- pmin(rep(coef_val, p), 5)
  names(coefficient_initial) <- colnames(cell_features)
  coefficient_initial["Intercept"] <- 0
  
  # apply loocv
  tmp <- sapply(1:loocv_len, function(i){
    print(paste0("Dropping lineage #", i, " out of ", loocv_len, " (", loocv_lineages[i], ")"))
    
    cell_features_train <- cell_features[-loocv_idx_list[[i]],var_names,drop = F]
    cell_lineage_train <- cell_lineage[-loocv_idx_list[[i]]]
    lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) == loocv_lineages[i])]
    
    set.seed(10)
    tmp_res <- multiomeFate:::lineage_imputation(cell_features = cell_features_train,
                                                 cell_lineage = cell_lineage_train,
                                                 coefficient_initial = coefficient_initial,
                                                 lineage_future_count = lineage_future_count_train,
                                                 random_initializations = 10,
                                                 verbose = 0)
    
    cell_features_test <- cell_features[loocv_idx_list[[i]],var_names,drop = F]
    cell_lineage_test <- cell_lineage[loocv_idx_list[[i]]]
    lineage_future_count_test <- lineage_future_count[loocv_lineages[i]]
    
    loglik_val <- multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_test,
                                                        cell_lineage = cell_lineage_test,
                                                        coefficient_vec = tmp_res$fit$coefficient_vec,
                                                        lineage_future_count = lineage_future_count_test)
    loglik_val
  })
  
  loocv_mat[,iteration] <- tmp
  
  save(coefficient_list_list, loocv_mat, date_of_run, session_info,
       file = paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_day0_lineage-imputation_stepdown_step2-tmp.RData"))
} 

save(date_of_run, session_info,
     cell_features_full,
     cell_lineage,
     coefficient_list_list, 
     lineage_current_count,
     lineage_future_count,
     loocv_mat,
     tab_mat,
     file = paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_day0_lineage-imputation_stepdown_step2.RData"))

print("Done! :)")