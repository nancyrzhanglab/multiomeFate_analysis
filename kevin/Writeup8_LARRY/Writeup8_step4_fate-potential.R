rm(list=ls())
library(Seurat)
library(multiomeFate)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step2_lineage-plotting.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# we'll focus on predict all the day6 cells at multiple timepoints from day4

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

for(treatment in treatment_vec){
  day_later <- treatment
  seurat_object2 <- seurat_object
  
  # keep only the relevant cells
  keep_vec <- rep(FALSE, ncol(seurat_object2))
  idx <- which(seurat_object2$time_celltype %in% day_early_vec)
  keep_vec[idx] <- TRUE
  seurat_object2$keep <- keep_vec
  seurat_object2 <- subset(seurat_object2, keep == TRUE)
  
  tab_mat <- table(seurat_object2$assigned_lineage, seurat_object2$time_celltype)
  
  # keep only the relevant cells for this analysis
  keep_vec <- rep(FALSE, ncol(seurat_object2))
  keep_vec[which(seurat_object2$time_celltype %in% day_early_vec)] <- TRUE
  seurat_object2$keep <- keep_vec
  seurat_object2 <- subset(seurat_object2, keep == TRUE)
  
  # construct cell_features matrix
  topic_mat <- seurat_object2[[paste0("fasttopic_", treatment)]]@cell.embeddings
  cell_features <- scale(topic_mat)
  p <- ncol(cell_features)
  
  cell_lineage <- seurat_object2$assigned_lineage
  uniq_lineage <- sort(unique(cell_lineage))
  lineage_current_count <- tab_mat[uniq_lineage, day_early_vec]
  lineage_future_count <- tab_mat[uniq_lineage, day_later]
  tab_mat <- tab_mat[uniq_lineage,]
  
  lambda_initial <- 3
  
  ###############################
  
  # construct the folds
  num_folds <- 10
  lineages_ordered <- rownames(tab_mat)[order(tab_mat[,day_later], decreasing = T)]
  num_lineages <- length(lineages_ordered)
  num_per_fold <- ceiling(length(lineages_ordered)/num_folds)
  set.seed(10)
  for(i in 1:num_per_fold){
    idx_vec <- ((i-1)*num_per_fold+1):min(i*num_per_fold, num_lineages)
    lineages_ordered[sample(idx_vec)] <- lineages_ordered[idx_vec]
  }
  fold_lineage_list <- lapply(1:num_folds, function(i){
    lineages_ordered[unique(pmin(i + (0:num_per_fold)*num_folds, num_lineages))]
  })
  names(fold_lineage_list) <- paste0("fold:", 1:num_folds)
  
  cv_cell_list <- lapply(fold_lineage_list, function(lineages){
    unlist(lapply(lineages, function(lineage){
      which(cell_lineage == lineage)
    }))
  })
  names(cv_cell_list) <- names(fold_lineage_list)
  
  # start the fitting process
  
  cv_fit_list <- vector("list", length = num_folds)
  names(cv_fit_list) <- names(fold_lineage_list)
  
  set.seed(10)
  for(i in 1:num_folds){
    fold <- names(fold_lineage_list)[i]
    print(paste0("Dropping fold #", i, " out of ", num_folds))
    
    #################
    
    # training 
    cell_features_train <- cell_features[-cv_cell_list[[fold]],,drop = F]
    cell_lineage_train <- cell_lineage[-cv_cell_list[[fold]]]
    lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) %in% fold_lineage_list[[fold]])]
    
    set.seed(10)
    train_fit <- multiomeFate:::lineage_imputation_sequence(
      cell_features = cell_features_train,
      cell_lineage = cell_lineage_train,
      lambda_initial = lambda_initial,
      lambda_sequence_length = 50,
      lineage_future_count = lineage_future_count_train,
      verbose = 1
    )
    
    #################
    
    # training evaluation
    lambda_sequence <- train_fit$lambda_sequence
    
    train_loglik <- sapply(1:length(lambda_sequence), function(kk){
      multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_train,
                                            cell_lineage = cell_lineage_train,
                                            coefficient_vec = train_fit$fit_list[[kk]]$coefficient_vec,
                                            lineage_future_count = lineage_future_count_train,
                                            lambda = 0)
    })
    
    #################
    
    # testing
    cell_features_test <- cell_features[cv_cell_list[[fold]],,drop = F]
    cell_lineage_test <- cell_lineage[cv_cell_list[[fold]]]
    lineage_future_count_test <- lineage_future_count[which(names(lineage_future_count) %in% fold_lineage_list[[fold]])]
    
    test_loglik <- sapply(1:length(lambda_sequence), function(kk){
      multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_test,
                                            cell_lineage = cell_lineage_test,
                                            coefficient_vec = train_fit$fit_list[[kk]]$coefficient_vec,
                                            lineage_future_count = lineage_future_count_test,
                                            lambda = 0)
    })
    
    #################
    
    cv_fit_list[[fold]] <- list(test_loglik = test_loglik,
                                train_loglik = train_loglik,
                                train_fit = train_fit)
    
    save(cv_fit_list, date_of_run, session_info,
         file = paste0("../../../../out/kevin/Writeup6r/Writeup6r_", treatment, "_", day_early, "_lineage-imputation-tmp.RData"))
  }
  
  
  print(paste0("Finished ", treatment))
  date_of_run <- Sys.time()
  session_info <- devtools::session_info()
  save(date_of_run, session_info,
       cell_features,
       cell_lineage,
       cv_cell_list,
       cv_fit_list, 
       fold_lineage_list,
       lineage_current_count,
       lineage_future_count,
       tab_mat,
       treatment,
       file = paste0("~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step5_growth-potential_", treatment, ".RData"))
  
  print("=========")
}

print("Done! :)")
