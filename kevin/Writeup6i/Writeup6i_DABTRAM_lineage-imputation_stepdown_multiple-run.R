# try data fission
rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_features <- cbind(1, scale(all_data$pseudotime[colnames(all_data2)]))
colnames(cell_features) <- c("Intercept", "Pseudotime")
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]
lineage_future_count <- tab_mat[uniq_lineage,"week5_DABTRAM"]

topic_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[names(cell_lineage),]
pseudotime_vec <- all_data$pseudotime[names(cell_lineage)]
atac_mat <- all_data2[["lsi"]]@cell.embeddings

log10pval_vec <- sapply(1:ncol(atac_mat), function(j){
  x <- atac_mat[names(cell_lineage),j]
  y <- as.factor(all_data2$tier_vec[names(cell_lineage)])
  res <- stats::oneway.test(x ~ y)
  -log10(res$p.value)
})
names(log10pval_vec) <- colnames(atac_mat)

# let's try using all the topics #Clueless
cell_features <- cbind(1, scale(topic_mat), 
                       scale(atac_mat[names(cell_lineage),which(log10pval_vec>=18)]), 
                       scale(pseudotime_vec))
p <- ncol(cell_features)
colnames(cell_features)[1] <- "Intercept"
colnames(cell_features)[p] <- "Pseudotime"

tmp <- quantile(abs(cell_features[,-1]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
coefficient_initial <- c(0, rep(coef_val, p-1))
names(coefficient_initial) <- colnames(cell_features)

##################

# try to fission the cell_features
# https://arxiv.org/abs/2112.11079

n <- nrow(cell_features)
trials <- 10
noise_list <- lapply(1:trials, function(trial){
  set.seed(10*trial)
  
  z_mat <- matrix(0, nrow = n, ncol = ncol(cell_features))
  for(j in 2:ncol(z_mat)){
    sd_val <- stats::sd(cell_features[,j])
    z_mat[,j] <- stats::rnorm(n = n, mean = 0, sd = sd_val)
  }
  colnames(z_mat) <- colnames(cell_features)
  
  z_mat
})

coefficient_list_list <- vector("list", length = 1)
iteration <- 1
while(TRUE){
  print(paste0("On iteration ", iteration))
  
  if(iteration > 1) {
    var_keep <- coefficient_list[[iteration-1]]$variables_to_be_kept
  } else {
    var_keep <- colnames(cell_features)
  }
  print(var_keep)
  if(length(var_keep) < 2) break()
  
  p <- length(var_keep)
  tmp <- quantile(abs(cell_features[,setdiff(var_keep, "Intercept")]), probs = 0.95)
  coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
  coefficient_initial <- c(0, rep(coef_val, p-1))
  names(coefficient_initial) <- colnames(var_keep)
  
  # apply the trials
  coefficient_list <- lapply(1:trials, function(trial){
    print(paste0("Working on trial ", trial))
    training_mat <- cell_features[,var_keep] + noise_list[[trial]][,var_keep]
    testing_mat <- cell_features[,var_keep] - noise_list[[trial]][,var_keep]
    
    set.seed(10)
    lineage_res <- multiomeFate:::lineage_imputation(cell_features = training_mat,
                                                     cell_lineage = cell_lineage,
                                                     coefficient_initial = coefficient_initial,
                                                     lineage_future_count = lineage_future_count,
                                                     random_initializations = 10,
                                                     verbose = 1)
    coefficient_vec <- lineage_res$fit$coefficient_vec
   
    loglik_val <- multiomeFate:::evaluate_loglikelihood(cell_features = testing_mat,
                                                        cell_lineage = cell_lineage,
                                                        coefficient_vec = lineage_res$fit$coefficient_vec,
                                                        lineage_future_count = lineage_future_count)
    list(coefficient_vec = coefficient_vec,
         testing_obj = loglik_val,
         training_obj = lineage_res$fit$objective_val)
  })
  names(coefficient_list) <- paste0("trial_", 1:length(coefficient_list))
  
  coef_mat <- sapply(coefficient_list, function(x){abs(x$coefficient_vec)})
  coef_mean <- rowMeans(coef_mat)
  variable_to_remove <- rownames(coef_mat)[which.min(coef_mean)]
  variables_to_be_kept <- setdiff(rownames(coef_mat), variable_to_remove)
  coefficient_list$variables_to_be_kept <- variables_to_be_kept
  
  coefficient_list_list[[iteration]] <- coefficient_list
  iteration <- iteration+1
  
  save(coefficient_list_list, date_of_run, session_info,
       file = "../../../../out/kevin/Writeup6i/Writeup6i_DABTRAM_lineage-imputation_stepdown_multiple-run.RData")
}

save(coefficient_list_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6i/Writeup6i_DABTRAM_lineage-imputation_stepdown_multiple-run.RData")
