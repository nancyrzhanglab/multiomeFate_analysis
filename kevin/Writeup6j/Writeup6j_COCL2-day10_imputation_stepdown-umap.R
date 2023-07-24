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
var_names2 <- setdiff(var_names, "Intercept")

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

#####################

zz <- cell_features[,var_names2]
colnames(zz) <- paste0("dimred_", 1:ncol(zz))
all_data2[["dimred"]] <- Seurat::CreateDimReducObject(embeddings = zz,
                                                      key = "dimred")
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, 
                            reduction = "dimred",
                            reduction.name = "dimredUMAP",
                            dims = 1:19)
all_data2$imputed_count <- log(cell_imputed_count)

# grDevices::colorRampPalette("lightgray", "blue")
p1 <- scCustomize::FeaturePlot_scCustom(all_data2, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        na_cutoff = quantile(all_data2$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "dimredUMAP", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("COCL2 Day10 imputed counts\n(Stepup from RNA fasttopics, ATAC LSI), (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_COCL2-day10_imputation_stepdown_umap-onlyDay10.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

##################


cell_features <- cbind(scale(fasttopic_mat), 
                       scale(lsi_mat))
zz <- cell_features
colnames(zz) <- paste0("dimred_", 1:ncol(zz))
all_data2[["dimred"]] <- Seurat::CreateDimReducObject(embeddings = zz,
                                                      key = "dimred")
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, 
                             reduction = "dimred",
                             reduction.name = "dimredUMAP",
                             dims = 1:ncol(zz))
all_data2$imputed_count <- log(cell_imputed_count)

# grDevices::colorRampPalette("lightgray", "blue")
p1 <- scCustomize::FeaturePlot_scCustom(all_data2, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        na_cutoff = quantile(all_data2$imputed_count, probs = 0.05, na.rm = T),
                                        na_color = "bisque",
                                        reduction = "dimredUMAP", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("COCL2 Day10 imputed counts\n(All RNA fasttopics, ATAC LSI), (Log-scale)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_COCL2-day10_imputation_stepdown_umap-onlyDay10_allvars.png"),
                p1, device = "png", width = 7, height = 5, units = "in")


