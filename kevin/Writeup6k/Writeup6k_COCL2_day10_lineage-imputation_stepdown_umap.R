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

cell_imputed_count <- as.numeric(exp(cell_features %*% lineage_res$coefficient_vec))
names(cell_imputed_count) <- rownames(cell_features)

lineage_imputed_count <- sapply(uniq_lineage, function(lineage){
  sum(cell_imputed_count[which(cell_lineage == lineage)])
})

imputed_vec <- rep(NA, ncol(all_data))
names(imputed_vec) <- colnames(all_data)
imputed_vec[names(cell_imputed_count)] <- log10(cell_imputed_count)
all_data$imputed_count <- imputed_vec

###################

keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset == "week5_COCL2")] <- TRUE
keep_vec[intersect(
  which(all_data$dataset == "day10_COCL2"),
  which(all_data$imputed_count >= -1)
)] <- TRUE
all_data$keep <- keep_vec

table(all_data$dataset, all_data$keep)
all_data2 <- subset(all_data, keep == TRUE)

all_data2[["pca"]] <- NULL
all_data2[["umap"]] <- NULL
Seurat::DefaultAssay(all_data2) <- "Saver"
set.seed(10)
all_data2 <- Seurat::RunPCA(all_data2, 
                           features = Seurat::VariableFeatures(object = all_data),
                           verbose = F)
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, 
                            reduction = "pca",
                            dims = 1:30)

# let's first plot the data as-is
col_vec <- c("darkblue", "#DE2D26")
names(col_vec) <- c("week5_COCL2", "day10_COCL2")
plot1 <- Seurat::DimPlot(all_data2, reduction = "umap", cols = col_vec,
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2 (High growth Day10 + Week5)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_day10highgrowth-week5_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# grDevices::colorRampPalette("lightgray", "blue")
p1 <- scCustomize::FeaturePlot_scCustom(all_data2, 
                                        colors_use = list("lightpink", "blue", "darkblue"),
                                        na_cutoff = NULL,
                                        na_color = "bisque",
                                        reduction = "umap", 
                                        features = "imputed_count")
p1 <- p1 + ggplot2::ggtitle("COCL2 Day10 imputed counts\n(Only high growth Day10's)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_day10highgrowth-week5_umap_day10-highlight.png"),
                p1, device = "png", width = 7, height = 5, units = "in")




