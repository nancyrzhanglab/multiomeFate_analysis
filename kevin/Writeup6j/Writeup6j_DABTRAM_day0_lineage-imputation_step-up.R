rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM-day0_extracted.RData")
all_data2$tier_vec <- all_data2$keep

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]

# construct cell_features matrix
topic_mat <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[names(cell_lineage),]
atac_mat <- all_data2[["lsi"]]@cell.embeddings

log10pval_vec <- sapply(1:ncol(atac_mat), function(j){
  x <- atac_mat[names(cell_lineage),j]
  y <- as.factor(all_data2$tier_vec[names(cell_lineage)])
  res <- stats::oneway.test(x ~ y)
  -log10(res$p.value)
})
names(log10pval_vec) <- colnames(atac_mat)

# let's try using all the topics #Clueless
cell_features_full <- cbind(1, scale(topic_mat), 
                            scale(atac_mat[names(cell_lineage),which(log10pval_vec>=2)]))
p <- ncol(cell_features_full)
colnames(cell_features_full)[1] <- "Intercept"

cell_features_full <- cell_features_full[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day0"]
lineage_future_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]

tmp <- quantile(abs(cell_features_full[,-1]), probs = 0.95)
coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
coefficient_initial <- rep(coef_val, p)
names(coefficient_initial) <- colnames(cell_features_full)
coefficient_initial["Intercept"] <- 0

#######################

coefficient_list_list <- vector("list", length = 1)
var_current <- "Intercept"
iteration <- 1
while(TRUE){
  print(paste0("On iteration ", iteration))
  
  if(iteration > 1) var_current <- coefficient_list_list[[iteration-1]]$var_next_iteration
  var_try <- setdiff(colnames(cell_features_full), var_current)
  if(length(var_try) <= 1) break()
  
  attempt_vec <- sapply(var_try, function(variable){
    print(variable)
    cell_features <- cell_features_full[,c(var_current, variable),drop=F]
    p <- ncol(cell_features)
    tmp <- quantile(abs(cell_features[,-which(colnames(cell_features) == "Intercept"),drop=F]), probs = 0.95)
    coef_val <- 2*log(max(lineage_future_count/lineage_current_count))/((p-1)*tmp)
    coefficient_initial <- rep(coef_val, p)
    names(coefficient_initial) <- colnames(cell_features)
    coefficient_initial["Intercept"] <- 0
    
    set.seed(10)
    tmp_res <- multiomeFate:::lineage_imputation(cell_features = cell_features,
                                                 cell_lineage = cell_lineage,
                                                 coefficient_initial = coefficient_initial,
                                                 lineage_future_count = lineage_future_count,
                                                 random_initializations = 10,
                                                 verbose = 0)
    
    tmp_res$fit$objective_val
  })
  names(attempt_vec) <- var_try
  
  # pick the variable to add based on the best fit (i.e., smallest objective)
  var_selected <- names(attempt_vec)[which.min(attempt_vec)]
  
  coefficient_list_list[[iteration]] <- list(
    attempt_vec = attempt_vec,
    var_next_iteration = sort(c(var_current, var_selected)),
    var_selected = var_selected
  )
  iteration <- iteration+1
  
  save(coefficient_list_list, date_of_run, session_info,
       file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_step-up.RData")
}

save(coefficient_list_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_step-up.RData")

