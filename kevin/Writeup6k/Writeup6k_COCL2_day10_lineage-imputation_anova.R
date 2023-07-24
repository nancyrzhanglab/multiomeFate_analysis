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

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_COCL2",
                            dims = 1:30)

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

#################

tab_mat <- tab_mat[order(tab_mat[,"week5_COCL2"], decreasing = T),]
lineage_idx <- intersect(which(tab_mat[,"day10_COCL2"] >= 3),
                         which(tab_mat[,"week5_COCL2"] >= 10))
large_lineages <- rownames(tab_mat)[lineage_idx]
all_data3 <- subset(all_data, 
                    dataset == "day10_COCL2" & assigned_lineage %in% large_lineages)

large_lineages2 <- large_lineages
lineage_vec <- all_data3$assigned_lineage

######################

df <- data.frame(lineage = lineage_vec,
                 imputed_count = all_data3$imputed_count)
df$lineage <- as.factor(df$lineage)

set.seed(10)
anova_res <- stats::oneway.test(imputed_count ~ lineage, data = df)

