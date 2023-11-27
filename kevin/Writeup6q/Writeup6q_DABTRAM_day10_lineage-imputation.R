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
treatment <- "DABTRAM"

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

# keep only the relevant cells for this analysis
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset == paste0("day10_", treatment))] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

# construct cell_features matrix
topic_mat <- all_data[[paste0("fasttopic_", treatment)]]@cell.embeddings
atac_mat <- all_data[[paste0("peakVI_", treatment)]]@cell.embeddings

# let's try using all the topics
cell_features_full <- cbind(1, scale(topic_mat), scale(atac_mat))
p <- ncol(cell_features_full)
colnames(cell_features_full)[1] <- "Intercept"

cell_lineage <- all_data$assigned_lineage
uniq_lineage <- sort(unique(cell_lineage))
lineage_current_count <- tab_mat[uniq_lineage, paste0("day10_", treatment)]
lineage_future_count <- tab_mat[uniq_lineage, paste0("week5_", treatment)]

tmp_res <- .compute_initial_parameters(cell_features = cell_features,
                                       cell_lineage = cell_lineage,
                                       lineage_future_count = lineage_future_count,
                                       multipler = 1e4)
lambda_initial <- tmp_res$lambda_initial

loocv_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 20),
                                              which(tab_mat[,paste0("day10_", treatment)] >= 1))]
loocv_len <- length(loocv_lineages)
loocv_idx_list <- lapply(loocv_lineages, function(lineage){
  which(cell_lineage == lineage)
})

#########################

loocv_fit_list <- numeric("list", length = loocv_len)
names(loocv_fit_list) <- loocv_lineages

for(i in 1:loocv_len){
  print(paste0("Dropping lineage #", i, " out of ", loocv_len, " (", loocv_lineages[i], ")"))
  
  cell_features_train <- cell_features[-loocv_idx_list[[i]],var_names,drop = F]
  cell_lineage_train <- cell_lineage[-loocv_idx_list[[i]]]
  lineage_future_count_train <- lineage_future_count[-which(names(lineage_future_count) == loocv_lineages[i])]
  
  set.seed(10)
  loocv_fit_list[[loocv_lineages[i]]] <- multiomeFate:::lineage_imputation_sequence(
    cell_features = cell_features_train,
    cell_lineage = cell_lineage_train,
    lambda_initial = lambda_initial,
    lineage_future_count = lineage_future_count_train,
    verbose = 0
  )
  
  save(loocv_fit_list, date_of_run, session_info,
       file = paste0("../../../../out/kevin/Writeup6q/Writeup6q_", treatment, "_day10_lineage-imputation-tmp.RData"))
}

save(date_of_run, session_info,
     cell_features_full,
     cell_lineage,
     lineage_current_count,
     lineage_future_count,
     loocv_fit_list, 
     tab_mat,
     file = paste0("../../../../out/kevin/Writeup6q/Writeup6q_", treatment, "_day10_lineage-imputation.RData"))

print("Done! :)")
