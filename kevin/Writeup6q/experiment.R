rm(list=ls())
library(Seurat)
library(multiomeFate)

load("../../../../out/kevin/Writeup6p/Writeup6p_all-data_lightweight_noATAC.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)
treatment <- "CIS"

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

# keep only the relevant cells for this analysis
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset == "day0")] <- TRUE
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
lineage_current_count <- tab_mat[uniq_lineage, "day0"]
lineage_future_count <- tab_mat[uniq_lineage, paste0("day10_", treatment)]

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
  lineage_future_count = lineage_future_count,
  multipler = 1e4
)
lambda_initial <- tmp_res$lambda_initial

loocv_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("day10_", treatment)] >= 20),
                                              which(tab_mat[,"day0"] >= 1))]
loocv_len <- length(loocv_lineages)
loocv_idx_list <- lapply(loocv_lineages, function(lineage){
  which(cell_lineage == lineage)
})
names(loocv_idx_list) <- loocv_lineages

cell_features_safe <- cell_features
cell_lineage_safe <- cell_lineage
lineage_future_count_safe <- lineage_future_count

cell_features <- cell_features_safe
cell_lineage <- cell_lineage_safe
lineage_future_count <- lineage_future_count_safe

#################

loocv_fit_list <- vector("list", length = loocv_len)
names(loocv_fit_list) <- loocv_lineages

res <- multiomeFate:::.compute_initial_parameters(cell_features = cell_features,
                                                  cell_lineage = cell_lineage,
                                                  lineage_future_count = lineage_future_count,
                                                  multipler = 100)

i <- 1
lineage <- loocv_lineages[i]

cell_features_train <- cell_features[-loocv_idx_list[[lineage]],,drop = F]
cell_lineage_train <- cell_lineage[-loocv_idx_list[[lineage]]]
lineage_future_count_train <- lineage_future_count[which(names(lineage_future_count) != lineage)]

cell_features_test <- cell_features[loocv_idx_list[[lineage]],,drop = F]
cell_lineage_test <- cell_lineage[loocv_idx_list[[lineage]]]
lineage_future_count_test <- lineage_future_count[lineage]

set.seed(10)
lambda <- 10000
imputation_res1 <- multiomeFate:::lineage_imputation(cell_features = cell_features_train,
                                                     cell_lineage = cell_lineage_train,
                                                     coefficient_initial_list = res$coefficient_initial,
                                                     lineage_future_count = lineage_future_count_train,
                                                     lambda = lambda,
                                                     random_initializations = 10,
                                                     upper_randomness = 5,
                                                     verbose = 0)

set.seed(10)
lambda <- 10
imputation_res2 <- multiomeFate:::lineage_imputation(cell_features = cell_features_train,
                                                     cell_lineage = cell_lineage_train,
                                                     coefficient_initial_list = res$coefficient_initial,
                                                     lineage_future_count = lineage_future_count_train,
                                                     lambda = lambda,
                                                     random_initializations = 10,
                                                     upper_randomness = 5,
                                                     verbose = 0)

set.seed(10)
lambda <- 0
imputation_res3 <- multiomeFate:::lineage_imputation(cell_features = cell_features_train,
                                                     cell_lineage = cell_lineage_train,
                                                     coefficient_initial_list = res$coefficient_initial,
                                                     lineage_future_count = lineage_future_count_train,
                                                     lambda = lambda,
                                                     random_initializations = 10,
                                                     upper_randomness = 5,
                                                     verbose = 0)

round(cbind(imputation_res1$fit$coefficient_vec,
            imputation_res2$fit$coefficient_vec,
            imputation_res3$fit$coefficient_vec), 2)

coef_list <- list(imputation_res1$fit$coefficient_vec,
                  imputation_res2$fit$coefficient_vec,
                  imputation_res3$fit$coefficient_vec)

########################

sapply(1:3, function(j){
  multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_test,
                                        cell_lineage = cell_lineage_test,
                                        coefficient_vec = coef_list[[j]],
                                        lineage_future_count = lineage_future_count_test,
                                        lambda = lambda)
})

sapply(1:3, function(j){
  multiomeFate:::evaluate_loglikelihood(cell_features = cell_features_train,
                                        cell_lineage = cell_lineage_train,
                                        coefficient_vec = coef_list[[j]],
                                        lineage_future_count = lineage_future_count_train,
                                        lambda = lambda)
})

########################

tmp <- multiomeFate:::.lineage_cleanup(cell_features = cell_features_train,
                                       cell_lineage = cell_lineage_train,
                                       lineage_future_count = lineage_future_count_train)

multiomeFate:::.lineage_objective(tmp$cell_features,
                                  tmp$cell_lineage,
                                  tmp$cell_lineage_idx_list,
                                  imputation_res1$fit$coefficient_vec,
                                  lambda,
                                  tmp$lineage_future_count)

######################

cell_features <- tmp$cell_features
cell_lineage <- tmp$cell_lineage
cell_lineage_idx_list <- tmp$cell_lineage_idx_list
lineage_future_count <- tmp$lineage_future_count
coefficient_vec <- imputation_res1$fit$coefficient_vec

uniq_lineages <- names(lineage_future_count)
cell_names <- rownames(cell_features)

scalar1 <- as.numeric(exp(cell_features %*% coefficient_vec))
names(scalar1) <- cell_names
scalar2 <- sapply(uniq_lineages, function(lineage){
  log(sum(scalar1[cell_lineage_idx_list[[lineage]]]))
})

idx_notintercept <- which(names(coefficient_vec) != "Intercept")
scalar3 <- multiomeFate:::.l2norm(coefficient_vec[idx_notintercept])

sum(scalar1) - sum(lineage_future_count*scalar2) + lambda*scalar3^2

