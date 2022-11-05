rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)

load("../../../../out/kevin/Writeup5a/Writeup5a_tcca_RNA-geneActivity.RData")
source("../Writeup5a/color_palette.R")

#CIS: 3,4,7,8,14,18,19,21,23,27,27,29
#COCL2: 1,2,4,5,7,11,12,13,14,15,19,21,22,26,27,28,30
#DABTRAM: 1,5,7,8,11,12,13,15,17,23,24,25,28,30

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
topic_list <- list("CIS" = c(3,4,7,8,14,18,19,21,23,27,27,29),
                   "COCL2" = c(1,2,4,5,7,11,12,13,14,15,19,21,22,26,27,28,30),
                   "DABTRAM" = c(1,5,7,8,11,12,13,15,17,23,24,25,28,30))
gene_all_list <- lapply(treatment_vec, function(treatment){
  topic_gene_mat <- all_data[[paste0("fasttopic_", treatment)]]@feature.loadings
  topic_vec <- topic_list[[treatment]]
  
  p <- nrow(topic_gene_mat)
  gene_list <- lapply(topic_vec, function(topic){
    rownames(topic_gene_mat)[which(topic_gene_mat[,topic] >= 15*1/p)]
  })
  gene_vec <- sort(unique(unlist(gene_list)))
  gene_vec
})
gene_all_vec <- sort(unique(unlist(gene_all_list)))
gene_all_vec <- intersect(intersect(gene_all_vec, rownames(multiSVD_obj$svd_1$v)), rownames(multiSVD_obj$svd_2$v))
length(gene_all_vec)

cell_cycling <- c(cc.genes$s.genes, cc.genes$g2m.genes)
length(intersect(cell_cycling, gene_all_vec))

#####################

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

######################

multiSVD_obj[["common_mat_1"]] <- NULL
multiSVD_obj[["distinct_mat_1"]] <- NULL
multiSVD_obj[["common_dimred_2"]] <- NULL
multiSVD_obj[["distinct_dimred_2"]] <- NULL

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

common_mat_1 <- multiSVD_obj$common_mat_1[,gene_all_vec]
common_mat_2 <- multiSVD_obj$common_mat_2[,gene_all_vec]

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
common_mat <- common_mat[-na_idx,]

mat_avg_split <- Matrix::t(sapply(1:length(combined_idx_list), function(i){
  if(i %% floor(length(combined_idx_list)/10) == 0) cat('*')
  matrixStats::colMedians(common_mat[combined_idx_list[[i]],,drop = F])
}))
rownames(mat_avg_split) <- names(combined_idx_list)

# mat_avg_split2 <- mat_avg_split[,selection_res$selected_variables]
set.seed(10)
svd_res <- irlba::irlba(mat_avg_split, nv = 20)
dimred <- svd_res$u %*% diag(svd_res$d)
rownames(dimred) <- rownames(mat_avg_split)

set.seed(10)
umap_res <- Seurat::RunUMAP(dimred)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- rownames(mat_avg_split)

lineage_vec2 <- sapply(rownames(mat_avg_split), function(x){
  strsplit(x, split = ":")[[1]][1]
})
dataset_vec2 <- sapply(rownames(mat_avg_split), function(x){
  strsplit(x, split = ":")[[1]][2]
})
col_vec <- sapply(dataset_vec2, function(x){
  col_palette[x]
})

set.seed(10)
idx <- sample(1:nrow(umap_mat))

png("../../../../out/figures/Writeup6b/Writeup6b_lineage-dataset_tcca_customGene2_UMAP_all.png",
    height = 3000, width = 3000, res = 500, units = "px")
par(mar = c(4,4,4,0.5))
plot(umap_mat[idx,1], umap_mat[idx,2],
     col = col_vec[idx], pch = 16, cex = 1, 
     xlab = "UMAP 1", ylab = "UMAP 2",
     main = "UMAP of lineage-dataset averages",
     xaxt = "n", yaxt = "n", bty = "n")
axis(1); axis(2)
graphics.off()

#################

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
treatment_vec <- c("CIS", "COCL2", "DABTRAM")
n <- nrow(umap_mat)
# day0
for(treatment in treatment_vec){
  lin_col <- rep(NA, n)
  names(lin_col) <- rownames(umap_mat)
  day10_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_",treatment)] >= 20)]
  week5_names <- rownames(tab_mat)[which(tab_mat[,paste0("week5_",treatment)] >= 20)]
  
  lin_col[intersect(names(lin_col), paste0(day10_names, ":day0"))] <- 3
  lin_col[intersect(names(lin_col), paste0(week5_names, ":day0"))] <- 2
  
  set.seed(10)
  idx <- sample(1:nrow(umap_mat))
  
  png(paste0("../../../../out/figures/Writeup6b/Writeup6b_lineage-dataset_tcca_customGene2_UMAP_day0_expandingDay10vsWeek5-", treatment, ".png"),
      height = 3000, width = 3000, res = 500, units = "px")
  par(mar = c(4,4,4,0.5))
  plot(umap_mat[idx,1], umap_mat[idx,2],
       col = col_vec[idx], pch = 16, cex = 1, 
       xlab = "UMAP 1", ylab = "UMAP 2",
       main = "UMAP of lineage-dataset averages",
       xaxt = "n", yaxt = "n", bty = "n")
  
  idx <- which(lin_col == 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = "white", pch = 16, cex = 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = lin_col[idx], pch = 16, cex = 2)
  
  idx <- which(lin_col == 2)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = "white", pch = 16, cex = 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = lin_col[idx], pch = 16, cex = 2)
  axis(1); axis(2)
  
  legend("bottomright", title=treatment,
         c("Day10>=20", "Week5>=20"), fill=c(3,2), cex = 0.8)
  graphics.off()
}

# day10
for(treatment in treatment_vec){
  lin_col <- rep(NA, n)
  names(lin_col) <- rownames(umap_mat)
  day10_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_",treatment)] >= 20)]
  week5_names <- rownames(tab_mat)[which(tab_mat[,paste0("week5_",treatment)] >= 20)]
  
  lin_col[intersect(names(lin_col), paste0(day10_names, ":day10_", treatment))] <- 3
  lin_col[intersect(names(lin_col), paste0(week5_names, ":day10_", treatment))] <- 2
  
  set.seed(10)
  idx <- sample(1:nrow(umap_mat))
  
  png(paste0("../../../../out/figures/Writeup6b/Writeup6b_lineage-dataset_tcca_customGene2_UMAP_day10_expandingDay10vsWeek5-", treatment, ".png"),
      height = 3000, width = 3000, res = 500, units = "px")
  par(mar = c(4,4,4,0.5))
  plot(umap_mat[idx,1], umap_mat[idx,2],
       col = col_vec[idx], pch = 16, cex = 1, 
       xlab = "UMAP 1", ylab = "UMAP 2",
       main = "UMAP of lineage-dataset averages",
       xaxt = "n", yaxt = "n", bty = "n")
  
  idx <- which(lin_col == 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = "white", pch = 16, cex = 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = lin_col[idx], pch = 16, cex = 2)
  
  idx <- which(lin_col == 2)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = "white", pch = 16, cex = 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = lin_col[idx], pch = 16, cex = 2)
  axis(1); axis(2)
  
  legend("bottomright", title=treatment,
         c("Day10>=20", "Week5>=20"), fill=c(3,2), cex = 0.8)
  graphics.off()
}

# week5
for(treatment in treatment_vec){
  lin_col <- rep(NA, n)
  names(lin_col) <- rownames(umap_mat)
  day10_names <- rownames(tab_mat)[which(tab_mat[,paste0("day10_",treatment)] >= 20)]
  week5_names <- rownames(tab_mat)[which(tab_mat[,paste0("week5_",treatment)] >= 20)]
  
  lin_col[intersect(names(lin_col), paste0(day10_names, ":week5_", treatment))] <- 3
  lin_col[intersect(names(lin_col), paste0(week5_names, ":week5_", treatment))] <- 2
  
  set.seed(10)
  idx <- sample(1:nrow(umap_mat))
  
  png(paste0("../../../../out/figures/Writeup6b/Writeup6b_lineage-dataset_tcca_customGene2_UMAP_week5_expandingDay10vsWeek5-", treatment, ".png"),
      height = 3000, width = 3000, res = 500, units = "px")
  par(mar = c(4,4,4,0.5))
  plot(umap_mat[idx,1], umap_mat[idx,2],
       col = col_vec[idx], pch = 16, cex = 1, 
       xlab = "UMAP 1", ylab = "UMAP 2",
       main = "UMAP of lineage-dataset averages",
       xaxt = "n", yaxt = "n", bty = "n")
  
  idx <- which(lin_col == 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = "white", pch = 16, cex = 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = lin_col[idx], pch = 16, cex = 2)
  
  idx <- which(lin_col == 2)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = "white", pch = 16, cex = 3)
  points(umap_mat[idx,1], umap_mat[idx,2],
         col = lin_col[idx], pch = 16, cex = 2)
  axis(1); axis(2)
  
  legend("bottomright", title=treatment,
         c("Day10>=20", "Week5>=20"), fill=c(3,2), cex = 0.8)
  graphics.off()
}


