# try data fission
rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_step-down_LOOCV.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis$fit$fit$coefficient_vec)
})
len <- length(coefficient_list_list)
train_vec <- sapply(coefficient_list_list, function(x){
  x$fit$fit$objective_val
})
test_vec <- sapply(1:len, function(i){
  if(i == 1) return(NA)
  removed_var <- setdiff(colnames(coefficient_list_list[[i]]$loocv_mat),
                         names(coefficient_list_list[[i]]$fit$fit$coefficient_vec))
  mean(coefficient_list_list[[i]]$loocv_mat[,removed_var])
})

####################################################

# let's try fitting the model with the selected model

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM-day0_extracted.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

idx <- which.min(test_vec)
lineage_res <- coefficient_list_list[[idx]]$fit
var_names <- names(lineage_res$fit$coefficient_vec)

fasttopic_mat <- all_data2[["fasttopic_DABTRAM"]]@cell.embeddings
lsi_mat <- all_data2[["lsi"]]@cell.embeddings

cell_features <- cbind(1, 
                       scale(fasttopic_mat), 
                       scale(lsi_mat))
colnames(cell_features)[1] <- "Intercept"
cell_features <- cell_features[,var_names]
cell_lineage <- all_data2$assigned_lineage
cell_lineage <- cell_lineage[!is.na(cell_lineage)]
cell_features <- cell_features[names(cell_lineage),]
uniq_lineage <- sort(unique(cell_lineage))
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_current_count <- tab_mat[uniq_lineage,"day0"]
lineage_future_count <- tab_mat[uniq_lineage,"day10_DABTRAM"]

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$fit$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)

lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log10(cell_imputed_count)

##################

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_DABTRAM",
                            dims = 1:30)

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log10(cell_imputed_count)

all_data$imputed_count <- imputed_vec

# grDevices::colorRampPalette("lightgray", "blue")
p1 <- scCustomize::FeaturePlot_scCustom(all_data, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        na_cutoff = quantile(all_data$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day0 imputed counts\n(Stepdown from RNA fasttopics, ATAC LSI), (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM-day0_imputation_stepdown-LOOCV_umap.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

####################
fit <- coefficient_list_list[[idx]]$fit
save(all_data, fit, lineage_future_count,
     date_of_run, session_info, 
     file = "../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day0_lineage-imputation_stepdown-LOOCV_postprocessed.RData")

