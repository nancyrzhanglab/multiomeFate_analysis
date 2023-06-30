rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2-day0_extracted.RData")
all_data2$tier_vec <- all_data2$keep

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]

# construct cell_features matrix
topic_mat <- all_data[["fasttopic_COCL2"]]@cell.embeddings[names(cell_lineage),]
atac_mat <- all_data2[["lsi"]]@cell.embeddings

log10pval_vec <- sapply(1:ncol(atac_mat), function(j){
  x <- atac_mat[names(cell_lineage),j]
  y <- as.factor(all_data2$tier_vec[names(cell_lineage)])
  res <- stats::oneway.test(x ~ y)
  -log10(res$p.value)
})
names(log10pval_vec) <- colnames(atac_mat)

# let's try using all the topics 
cell_features <- cbind(1, scale(topic_mat), 
                       scale(atac_mat[names(cell_lineage),which(log10pval_vec>=1)]))
p <- ncol(cell_features)
colnames(cell_features)[1] <- "Intercept"

cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day0"]
lineage_future_count <- tab_mat[uniq_lineage,"day10_COCL2"]

tmp <- quantile(abs(cell_features[,-1]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
coefficient_initial <- rep(coef_val, p)
names(coefficient_initial) <- colnames(cell_features)
coefficient_initial["Intercept"] <- 0

loocv_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,"day10_COCL2"] >= 25),
                                              which(tab_mat[,"day0"] >= 1))]
loocv_len <- length(loocv_lineages)
loocv_idx_list <- lapply(loocv_lineages, function(lineage){
  which(cell_lineage == lineage)
})

#########################

coefficient_list_list <- vector("list", length = 1)
iteration <- 1
while(TRUE){
  print(paste0("On iteration ", iteration))
  
  if(iteration > 1) {
    var_keep <- coefficient_list_list[[iteration-1]]$variables_to_be_kept
    var_rm <- coefficient_list_list[[iteration-1]]$variables_to_be_rm
  } else {
    var_keep <- colnames(cell_features)
    var_rm <- NULL
  }
  print(var_keep)
  if(length(var_keep) < 2) break()
  
  p <- length(var_keep)
  tmp <- quantile(abs(cell_features[,setdiff(var_keep, "Intercept")]), probs = 0.95)
  coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
  coefficient_initial <- rep(coef_val, p)
  names(coefficient_initial) <- var_keep
  coefficient_initial["Intercept"] <- 0
  
  if(iteration > 1){
    stopifnot(!is.null(var_rm), !is.na(var_rm), length(var_rm) >= 1)
    # for each variable being dropped:
    # compute their LOOCV among all the loocv_lineages
    loocv_mat <- matrix(NA, nrow = loocv_len, ncol = length(var_rm))
    rownames(loocv_mat) <- loocv_lineages
    colnames(loocv_mat) <- var_rm
    
    for(variable in var_rm){
      var_keep_tmp <- setdiff(var_keep, variable)
      
      # apply loocv
      tmp <- sapply(1:loocv_len, function(i){
        print(paste0("Dropping lineage #", i, " out of ", loocv_len, " (", loocv_lineages[i], ")"))
        
        cell_features_train <- cell_features[-loocv_idx_list[[i]],var_keep_tmp,drop = F]
        cell_lineage_train <- cell_lineage[-loocv_idx_list[[i]]]
        coefficient_initial_train <- coefficient_initial[var_keep_tmp]
        lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) == loocv_lineages[i])]
        
        set.seed(10)
        tmp_res <- multiomeFate:::lineage_imputation(cell_features = cell_features_train,
                                                     cell_lineage = cell_lineage_train,
                                                     coefficient_initial = coefficient_initial_train,
                                                     lineage_future_count = lineage_future_count_train,
                                                     random_initializations = 10,
                                                     verbose = 0)
        
        cell_features_test <- cell_features[loocv_idx_list[[i]],var_keep_tmp,drop = F]
        cell_lineage_test <- cell_lineage[loocv_idx_list[[i]]]
        lineage_future_count_test <- lineage_future_count[loocv_lineages[i]]
        
        loglik_val <- multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_test,
                                                            cell_lineage = cell_lineage_test,
                                                            coefficient_vec = tmp_res$fit$coefficient_vec,
                                                            lineage_future_count = lineage_future_count_test)
        loglik_val
      })
      
      loocv_mat[,variable] <- tmp
      
      save(coefficient_list_list, loocv_mat, date_of_run, session_info,
           file = "../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day0_lineage-imputation_step-down_LOOCV.RData")
    } 
    
    # now determine which variable to drop
    var_kick <- colnames(loocv_mat)[which.max(colMeans(loocv_mat))]
    stopifnot(var_kick %in% var_keep)
    var_keep <- setdiff(var_keep, var_kick)
  } else {
    loocv_mat <- NULL
  }
  
  # after deciding which variable to drop (or the first iteration), do a full fit
  #  to figure out the coefficient vector to report
  set.seed(10)
  lineage_res <- multiomeFate:::lineage_imputation(cell_features = cell_features[,var_keep,drop=F],
                                                   cell_lineage = cell_lineage,
                                                   coefficient_initial = coefficient_initial[var_keep],
                                                   lineage_future_count = lineage_future_count,
                                                   random_initializations = 10,
                                                   verbose = 1)
  coefficient_vec <- lineage_res$fit$coefficient_vec
  
  coefficient_list <- vector("list", length = 0)
  coefficient_list$fit <- lineage_res
  coefficient_list$loocv_mat <- loocv_mat
  coefficient_list$coefficient_vec <- coefficient_vec
  coefficient_list$variables_to_be_kept <- var_keep
  
  # schedule the 10 variables with the lowest coefficients to be dropped next iteration
  tmp <- names(coefficient_vec)[order(abs(coefficient_vec), decreasing = F)[1:min(5, length(coefficient_vec)-1)]]
  coefficient_list$variables_to_be_rm <- setdiff(tmp, "Intercept")
  
  coefficient_list_list[[iteration]] <- coefficient_list
  iteration <- iteration+1
  
  if(length(coefficient_list$variables_to_be_rm) == 1) break()
  
  save(coefficient_list_list, loocv_mat, date_of_run, session_info,
       file = "../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day0_lineage-imputation_step-down_LOOCV.RData")
}

save(coefficient_list_list, loocv_mat, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day0_lineage-imputation_step-down_LOOCV.RData")


