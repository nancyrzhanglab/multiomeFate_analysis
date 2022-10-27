rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)

load("../../../../out/kevin/Writeup6/Writeup6_tcca_selected-genes.RData")
load("../../../../out/kevin/Writeup5a/Writeup5a_tcca_RNA-geneActivity.RData")
all_data2 <- all_data
load("../../../../out/kevin/Writeup6/Writeup6_all-data_lineage-assigned.RData")
source("../Writeup5a/color_palette.R")

lineage_vec <- all_data$assigned_lineage
dataset_vec <- all_data$dataset

na_idx <- which(is.na(lineage_vec))
lineage_vec <- lineage_vec[-na_idx]
dataset_vec <- dataset_vec[-na_idx]

n <- length(lineage_vec)
combined_vec <- sapply(1:n, function(i){
  paste0(lineage_vec[i], ":", dataset_vec[i])
})
unique_combined <- sort(unique(combined_vec))
combined_idx_list <- sapply(unique_combined, function(combined_name){
  which(combined_vec == combined_name)
})
names(combined_idx_list) <- unique_combined
quantile(table(sapply(combined_idx_list, length)))

mat <- Matrix::t(all_data[["Saver"]]@scale.data)
mat_avg_split <- Matrix::t(sapply(1:length(combined_idx_list), function(i){
  if(i %% floor(length(combined_idx_list)/10) == 0) cat('*')
  matrixStats::colMedians(mat[combined_idx_list[[i]],,drop = F])
}))
rownames(mat_avg_split) <- names(combined_idx_list)

# mat_avg_split2 <- mat_avg_split[,selection_res$selected_variables]
set.seed(10)
svd_res <- irlba::irlba(mat_avg_split, nv = 20)
dimred <- svd_res$u %*% diag(svd_res$d)

set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
umap_mat <- umap_res@cell.embeddings

lineage_vec2 <- sapply(rownames(mat_avg_split), function(x){
  strsplit(x, split = ":")[[1]][1]
})
dataset_vec2 <- sapply(rownames(mat_avg_split), function(x){
  strsplit(x, split = ":")[[1]][2]
})
col_vec <- sapply(dataset_vec2, function(x){
  col_palette[x]
})

png("../../../../out/figures/Writeup6b/Writeup6b_lineage-dataset_UMAP_all.png",
    height = 3000, width = 3000, res = 500, units = "px")
par(mar = c(4,4,4,0.5))
plot(umap_mat[,1], umap_mat[,2],
     col = col_vec, pch = 16, cex = 1, 
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "UMAP of lineage-dataset averages",
     xaxt = "n", yaxt = "n", bty = "n")
axis(1); axis(2)
graphics.off()

########################################

multiSVD_obj[["common_mat_1"]] <- NULL
multiSVD_obj[["distinct_mat_1"]] <- NULL
multiSVD_obj[["common_dimred_2"]] <- NULL
multiSVD_obj[["distinct_dimred_2"]] <- NULL

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

common_mat_1 <- multiSVD_obj$common_mat_1[,selection_res$selected_variables]
common_mat_2 <- multiSVD_obj$common_mat_2[,selection_res$selected_variables]

common_mat_1 <- tiltedCCA:::.normalize_svd(input_obj = common_mat_1,
                                           averaging_mat = NULL,
                                           center = T,
                                           dims = 1:50,
                                           normalize_row = T,
                                           normalize_singular_value = T,
                                           recenter = F,
                                           rescale = F,
                                           scale = T,
                                           scale_max = NULL)
common_mat_2 <- tiltedCCA:::.normalize_svd(input_obj = common_mat_2,
                                           averaging_mat = NULL,
                                           center = T,
                                           dims = 2:50,
                                           normalize_row = T,
                                           normalize_singular_value = T,
                                           recenter = F,
                                           rescale = F,
                                           scale = T,
                                           scale_max = NULL)
common_mat <- cbind(common_mat_1, common_mat_2)
mat_avg_split2 <- common_mat[,selection_res$selected_variables]

mat <- Matrix::t(all_data[["Saver"]]@scale.data)
mat_avg_split <- Matrix::t(sapply(1:length(combined_idx_list), function(i){
  if(i %% floor(length(combined_idx_list)/10) == 0) cat('*')
  matrixStats::colMedians(mat[combined_idx_list[[i]],,drop = F])
}))
rownames(mat_avg_split) <- names(combined_idx_list)

# mat_avg_split2 <- mat_avg_split[,selection_res$selected_variables]
set.seed(10)
svd_res <- irlba::irlba(mat_avg_split, nv = 20)
dimred <- svd_res$u %*% diag(svd_res$d)

set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
umap_mat <- umap_res@cell.embeddings

lineage_vec2 <- sapply(rownames(mat_avg_split), function(x){
  strsplit(x, split = ":")[[1]][1]
})
dataset_vec2 <- sapply(rownames(mat_avg_split), function(x){
  strsplit(x, split = ":")[[1]][2]
})
col_vec <- sapply(dataset_vec2, function(x){
  col_palette[x]
})

png("../../../../out/figures/Writeup6b/Writeup6b_lineage-dataset_UMAP_all.png",
    height = 3000, width = 3000, res = 500, units = "px")
par(mar = c(4,4,4,0.5))
plot(umap_mat[,1], umap_mat[,2],
     col = col_vec, pch = 16, cex = 1, 
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "UMAP of lineage-dataset averages",
     xaxt = "n", yaxt = "n", bty = "n")
axis(1); axis(2)
graphics.off()

########################################

saver_umap <- all_data[["saverumap"]]@cell.embeddings
lineage_vec <- all_data$assigned_lineage
dataset_vec <- all_data$dataset

# make 7 plots, each showing how many cells are unlabeled
treatment_vec <- sort(unique(dataset_vec))
xlim <- range(saver_umap[,1]); ylim <- range(saver_umap[,2])
for(treatment in treatment_vec){
  print(treatment)
  
  png(paste0("../../../../out/figures/Writeup6b/Writeup6b_Saver-UMAP_", treatment, ".png"),
      height = 2500, width = 5000, res = 500, units = "px")
  par(mar = c(4,4,4,0.5), mfrow = c(1,2))
  idx1 <- which(dataset_vec == treatment)
  idx2 <- which(dataset_vec == treatment & !is.na(lineage_vec))
  
  plot(saver_umap[,1], saver_umap[,2],
       col = rgb(0.8,0.8,0.8), pch = 16, cex = 1, 
       xlab = "UMAP 1", ylab = "UMAP 2",
       main = paste0("Saver UMAP of ", treatment, " (All)"),
       xaxt = "n", yaxt = "n", bty = "n")
  points(saver_umap[idx1,1], saver_umap[idx1,2],
         col = "white", pch = 16, cex = 1)
  points(saver_umap[idx1,1], saver_umap[idx1,2],
         col = rgb(0.5,0.5,0.5,0.05), pch = 16, cex = 0.5)
  axis(1); axis(2)
  
  plot(saver_umap[,1], saver_umap[,2],
       col = rgb(0.8,0.8,0.8), pch = 16, cex = 1, 
       xlab = "UMAP 1", ylab = "UMAP 2",
       main = paste0("(Only lineage-assigned): ", round(100*length(idx2)/length(idx1)), "%"),
       xaxt = "n", yaxt = "n", bty = "n")
  points(saver_umap[idx2,1], saver_umap[idx2,2],
         col = "white", pch = 16, cex = 1)
  points(saver_umap[idx2,1], saver_umap[idx2,2],
         col = rgb(0.5,0.5,0.5,0.05), pch = 16, cex = 0.5)
  axis(1); axis(2)
  
  graphics.off()
}


