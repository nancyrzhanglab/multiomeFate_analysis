rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# How many cells need to be seen to count as "win"?
min_cell_count <- 50
earlierday_cond_str <- "day10_DABTRAM"
later_cond <- "week5_DABTRAM"

lintab <- as.data.frame.matrix(table(all_data$assigned_lineage, all_data$dataset)) # lineage by condition table

sel.lineages <- rownames(lintab)[which(lintab[,later_cond] >= min_cell_count)]  # lineages that survived in given condition
cells.1 <- which(all_data$assigned_lineage %in% sel.lineages 
                 & all_data$assigned_posterior > 0.5
                 & all_data$dataset == earlierday_cond_str)

# compute synchrony scores

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
all_data$alignment <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(all_data$alignment), max(all_data$alignment), length.out = num_color)
color_vec <- sapply(all_data$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})
names(alignment_vec) <- colnames(all_data)

######################

cell_idx <- intersect(intersect(which(all_data$dataset == "day10_DABTRAM"),
                                which(all_data$assigned_posterior > 0.5)),
                      which(!is.na(all_data$assigned_lineage)))
lineage_vec <- all_data$assigned_lineage[cell_idx]

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
future_num_vec <- sapply(lineage_vec, function(lineage){
  tab_mat[lineage,"week5_DABTRAM"]
})

# for a gene, compute the regression line regressing alignment onto future fitness across all the cells
target_genes <- colnames(rna_common)

order_vec <- order(future_num_vec, decreasing = F)
rna_mat <- rna_common[cell_idx[order_vec],target_genes]
atac_mat <- atac_pred[cell_idx[order_vec],target_genes]

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))
diff_mat <- (rna_mat2 - atac_mat2)/rna_mat2
colnames(diff_mat) <- colnames(rna_mat2)
diff_mat <- abs(diff_mat[,target_genes])
diff_mat <- pmin(diff_mat, 15)

p <- ncol(diff_mat); n <- nrow(diff_mat)
pred_diff_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = diff_mat[,j], x = log1p(future_num_vec[order_vec]))
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- log1p(future_num_vec[order_vec])
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_diff_mat) <- colnames(diff_mat)

rsquare_value <- sapply(1:ncol(pred_diff_mat), function(j){
  tmp_df <- data.frame(y = diff_mat[,j], x = log1p(future_num_vec[order_vec]))
  
  tmp_lm <- stats::lm(y ~ x, data = tmp_df)
  summary(tmp_lm)$r.squared
  # stats::coef(tmp_lm)["x"]
})

save(pred_diff_mat, rsquare_value,
     rna_mat, atac_mat, future_num_vec, alignment_vec,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6c/Writeup6c_RNA-ATAC_DABTRAM_synchrony-fitness_scan.RData")

