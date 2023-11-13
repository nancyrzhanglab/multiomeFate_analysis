rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6o/Writeup6o_all-data_peakvi.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)
treatment <- "COCL2"

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

# construct cell_features matrix
topic_mat <- all_data[[paste0("fasttopic_", treatment)]]@cell.embeddings
atac_mat <- all_data[[paste0("peakVI_", treatment)]]@cell.embeddings

# let's try using all the topics
cell_features_full <- cbind(1, scale(topic_mat), scale(atac_mat))
p <- ncol(cell_features_full)
colnames(cell_features_full)[1] <- "Intercept"

cell_lineage <- all_data2$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage, paste0("day10_", treatment)]
lineage_future_count <- tab_mat[uniq_lineage, paste0("week5_", treatment)]

tmp <- quantile(abs(cell_features_full[,-1]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
coefficient_initial <- rep(coef_val, p)
names(coefficient_initial) <- colnames(cell_features_full)
coefficient_initial["Intercept"] <- 0

#######################

coefficient_list_list <- vector("list", length = 1)
var_current <- colnames(cell_features_full)
iteration <- 1
while(TRUE){
  print(paste0("On iteration ", iteration))
  
  if(iteration > 1) var_current <- coefficient_list_list[[iteration-1]]$var_next_iteration
  if(length(var_current) <= 2) break()
  
  var_attempt <- setdiff(var_current, "Intercept")
  attempt_list <- lapply(var_attempt, function(variable){
    var_try <- setdiff(var_current, variable)
    cell_features <- cell_features_full[,var_try,drop=F]
    p <- ncol(cell_features)
    tmp <- quantile(abs(cell_features[,-which(colnames(cell_features) == "Intercept"),drop=F]), probs = 0.95)
    coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
    coefficient_initial <- pmin(rep(coef_val, p), 5)
    names(coefficient_initial) <- colnames(cell_features)
    coefficient_initial["Intercept"] <- 0
    
    set.seed(10)
    tmp_res <- multiomeFate:::lineage_imputation(cell_features = cell_features,
                                                 cell_lineage = cell_lineage,
                                                 coefficient_initial = coefficient_initial,
                                                 lineage_future_count = lineage_future_count,
                                                 random_initializations = 10,
                                                 verbose = 0)
    
    tmp_res$fit
  })
  names(attempt_list) <- var_attempt
  
  attempt_vec <- sapply(attempt_list, function(x){
    x$objective_val
  })
  
  # pick the variable to add based on the best fit (i.e., smallest objective)
  var_selected <- names(attempt_vec)[which.min(attempt_vec)]
  
  coefficient_list_list[[iteration]] <- list(
    attempt_vec = attempt_vec,
    fit = attempt_list[[which.min(attempt_vec)]],
    var_next_iteration = setdiff(var_current, var_selected),
    var_selected = var_selected
  )
  iteration <- iteration+1
  
  save(coefficient_list_list, date_of_run, session_info,
       file = paste0("../../../../out/kevin/Writeup6p/Writeup6p_", treatment, "_day10_lineage-imputation_stepdown-tmp.RData"))
}

save(date_of_run, session_info,
     cell_features_full,
     cell_lineage,
     coefficient_list_list, 
     lineage_current_count,
     lineage_future_count,
     tab_mat,
     file = paste0("../../../../out/kevin/Writeup6p/Writeup6p_", treatment, "_day10_lineage-imputation_stepdown.RData"))

print("Done! :)")
