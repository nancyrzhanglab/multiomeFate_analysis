rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day10_lineage-imputation_stepdown_step2.RData")
load("../../../../out/kevin/Writeup6i/Writeup6i_COCL2-day10_extracted.RData")
all_data2$tier_vec <- all_data2$keep

# extract the average test and train
len_vec <- sapply(coefficient_list_list, function(lis){
  length(lis$var_next_iteration)
})
len <- length(coefficient_list_list)
train_vec <- sapply(coefficient_list_list, function(x){
  x$fit$objective_val
})
test_vec <- colMeans(loocv_mat)

idx <- which.min(test_vec)
lineage_res <- coefficient_list_list[[idx]]$fit
var_names <- names(lineage_res$coefficient_vec)

fasttopic_mat <- all_data2[["fasttopic_COCL2"]]@cell.embeddings
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
lineage_current_count <- tab_mat[uniq_lineage,"day10_COCL2"]
lineage_future_count <- tab_mat[uniq_lineage,"week5_COCL2"]

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)
lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log10(cell_imputed_count)
all_data$imputed_count <- imputed_vec

############################

keep_vec <- !is.na(all_data$imputed_count)
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep == TRUE)
quantile(all_data2$imputed_count)

level_vec <- rep(NA, ncol(all_data2))
level_vec[which(all_data2$imputed_count >= quantile(all_data2$imputed_count, probs = 0.75))] <- "high"
level_vec[which(all_data2$imputed_count <= quantile(all_data2$imputed_count, probs = 0.25))] <- "low"
all_data2$level <- level_vec
Seurat::Idents(all_data2) <- "level"

# Seurat::DefaultAssay(all_data2) <- "Saver"
# set.seed(10)
# de_res <- Seurat::FindMarkers(
#   object = all_data2,
#   ident.1 = "high",
#   ident.2 = "low",
#   features = all_data2[["Saver"]]@var.features,
#   only.pos = FALSE,
#   min.pct = 0,
#   logfc.threshold = 0,
#   verbose = T
# )
# 
# head(de_res)
# length(which(de_res$p_val_adj <= 0.05))
# rownames(de_res)[de_res$p_val_adj <= 0.05]

# pvalue_list <- lapply(all_data2[["Saver"]]@var.features, function(gene){
#   x_vec <- all_data2[["Saver"]]@scale.data[gene, which(all_data2$level == "high")]
#   y_vec <- all_data2[["Saver"]]@scale.data[gene, which(all_data2$level == "low")]
#   
#   if(diff(range(x_vec)) <= 1e-6 || diff(range(y_vec)) <= 1e-6 ) return(NA)
#   
#   ttest_res <- stats::t.test(x = x_vec, y = y_vec)
#   list(stat = ttest_res$statistic, pvalue = ttest_res$p.value)
# })
# names(pvalue_list) <- all_data2[["Saver"]]@var.features
# pvalue_list <- pvalue_list[sapply(pvalue_list, function(x){!all(is.na(x))})]
# quantile(sapply(pvalue_list, function(x){x$stat}))

cor_vec <- sapply(all_data2[["Saver"]]@var.features, function(gene){
  x_vec <- all_data2[["Saver"]]@scale.data[gene,]
  y_vec <- all_data2$imputed_count
  
  stats::cor(x_vec, y_vec, method = "spearman")
})
quantile(cor_vec, na.rm = T)

