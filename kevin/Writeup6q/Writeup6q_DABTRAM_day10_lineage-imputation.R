rm(list=ls())
library(Seurat)

load("../../../../out/kevin/Writeup6p/Writeup6p_all-data_lightweight_noATAC.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment <- "DABTRAM"
day_early <- "day10"
day_later <- "week5"
day_early_full <- paste0(day_early, "_", treatment)
day_later_full <- paste0(day_later, "_", treatment)

#######################################

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

# keep only the relevant cells for this analysis
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset == day_early_full)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

# construct cell_features matrix
topic_mat <- all_data[[paste0("fasttopic_", treatment)]]@cell.embeddings
atac_mat <- all_data[[paste0("peakVI_", treatment)]]@cell.embeddings

# let's try using all the topics
cell_features <- cbind(scale(topic_mat), scale(atac_mat))
p <- ncol(cell_features)

cell_lineage <- all_data$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))
lineage_current_count <- tab_mat[uniq_lineage, day_early_full]
lineage_future_count <- tab_mat[uniq_lineage, day_later_full]

# do some initializations
tmp <- multiomeFate:::.lineage_cleanup(cell_features = cell_features,
                                       cell_lineage = cell_lineage,
                                       lineage_future_count = lineage_future_count)
cell_features <- tmp$cell_features
cell_lineage <- tmp$cell_lineage
lineage_future_count <- tmp$lineage_future_count

tmp_res <- multiomeFate:::.compute_initial_parameters(
  cell_features = cell_features,
  cell_lineage = cell_lineage,
  lambda_max = 1e4,
  lambda_min = 1e6,
  lineage_future_count = lineage_future_count,
  multipler = 1e4
)
lambda_initial <- tmp_res$lambda_initial

loocv_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,day_later_full] >= 20),
                                              which(tab_mat[,day_early_full] >= 1))]
loocv_len <- length(loocv_lineages)
loocv_idx_list <- lapply(loocv_lineages, function(lineage){
  which(cell_lineage == lineage)
})
names(loocv_idx_list) <- loocv_lineages

#########################

loocv_fit_list <- vector("list", length = loocv_len)
names(loocv_fit_list) <- loocv_lineages

for(i in 1:loocv_len){
  lineage <- loocv_lineages[i]
  print(paste0("Dropping lineage #", i, " out of ", loocv_len, " (", loocv_lineages[i], ")"))
  
  #################
  
  # training 
  cell_features_train <- cell_features[-loocv_idx_list[[lineage]],,drop = F]
  cell_lineage_train <- cell_lineage[-loocv_idx_list[[lineage]]]
  lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) == lineage)]
  
  set.seed(10)
  train_fit <- multiomeFate:::lineage_imputation_sequence(
    cell_features = cell_features_train,
    cell_lineage = cell_lineage_train,
    lambda_initial = lambda_initial,
    lambda_sequence_length = 30,
    lineage_future_count = lineage_future_count_train,
    verbose = 1
  )
  
  # train_fit$lambda_sequence
  
  #################
  
  # training evaluation
  lambda_sequence <- train_fit$lambda_sequence
  
  train_loglik <- sapply(1:length(lambda_sequence), function(kk){
    multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_train,
                                          cell_lineage = cell_lineage_train,
                                          coefficient_vec = train_fit$fit_list[[kk]]$coefficient_vec,
                                          lineage_future_count = lineage_future_count_train,
                                          lambda = lambda_sequence[kk])
  })
  
  # round(sapply(train_fit$fit_list, function(x){ x$coefficient_vec }),2)
  
  #################
  
  # testing
  cell_features_test <- cell_features[loocv_idx_list[[lineage]],,drop = F]
  cell_lineage_test <- cell_lineage[loocv_idx_list[[lineage]]]
  lineage_future_count_test <- lineage_future_count[lineage]
  
  test_loglik <- sapply(1:length(lambda_sequence), function(kk){
    multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_test,
                                          cell_lineage = cell_lineage_test,
                                          coefficient_vec = train_fit$fit_list[[kk]]$coefficient_vec,
                                          lineage_future_count = lineage_future_count_test,
                                          lambda = lambda_sequence[kk])
  })
  
  #################
  
  loocv_fit_list[[lineage]] <- list(test_loglik = test_loglik,
                                    train_loglik = train_loglik,
                                    train_fit = train_fit)
  
  save(loocv_fit_list, date_of_run, session_info,
       file = paste0("../../../../out/kevin/Writeup6q/Writeup6q_", treatment, "_day0_lineage-imputation-tmp.RData"))
}

save(date_of_run, session_info,
     cell_features,
     cell_lineage,
     lineage_current_count,
     lineage_future_count,
     loocv_fit_list, 
     tab_mat,
     file = paste0("../../../../out/kevin/Writeup6q/Writeup6q_", treatment, "_day0_lineage-imputation.RData"))

print("Done! :)")
