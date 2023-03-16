rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-chromAct_COCL2.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

cell_idx <- intersect(intersect(which(all_data$dataset == "day10_COCL2"),
                                which(all_data$assigned_posterior > 0.5)),
                      which(!is.na(all_data$assigned_lineage)))

lineage_vec <- all_data$assigned_lineage[cell_idx]
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
future_num_vec <- sapply(lineage_vec, function(lineage){
  tab_mat[lineage,"week5_COCL2"]
})
future_num_vec <- log1p(future_num_vec)

surviving_lineage <- rownames(tab_mat)[tab_mat[,"week5_COCL2"] >= 50]
cell_idx2 <- intersect(which(all_data$dataset == "week5_COCL2"),
                       which(all_data$assigned_lineage %in% surviving_lineage))

num_colors <- 20
color_palette <- grDevices::colorRampPalette(c("gray", "coral2"))(num_colors)
color_breaks <- seq(min(future_num_vec), max(future_num_vec), length.out = num_colors)
color_vec <- sapply(future_num_vec, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

color_vec2 <- rep(3, length(cell_idx2))

########

rna_common <- multiSVD_obj$common_mat_1
atac_pred <- multiSVD_obj$common_mat_2

colname_vec <- intersect(colnames(rna_common), colnames(atac_pred))
rna_common <- rna_common[,colname_vec]
atac_pred <- atac_pred[,colname_vec]

p <- ncol(atac_pred)
cor_vec <- sapply(1:p, function(j){
  stats::cor(rna_common[c(cell_idx, cell_idx2),j], atac_pred[c(cell_idx, cell_idx2),j])
})
quantile(cor_vec)
gene_idx <- order(cor_vec, decreasing = T)[1:25]

set.seed(10)
rng_idx <- sample(1:length(c(cell_idx, cell_idx2)))

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_COCL2_RNA-chromAct_leafplot_withweek5_genes_highest-correlation.png"),
    height = 3000, width = 3000, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(5,5))
for(j in gene_idx){
  plot(x = atac_pred[c(cell_idx, cell_idx2)[rng_idx],j],
       y = rna_common[c(cell_idx, cell_idx2)[rng_idx],j], 
       pch = 16, col = c(color_vec, color_vec2)[rng_idx],
       xlab = "Chromatin act. value", ylab = "RNA value", 
       main = paste0(colnames(rna_common)[j], "\nCorr: ", round(cor_vec[j],2)),
       cex.main = 0.7)
}
graphics.off()

####################################

label_vec <- sapply(lineage_vec, function(i){
  tab_mat[i,"week5_COCL2"] >= 20
})

classification_vec <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  dat_mat <- cbind(atac_pred[cell_idx,j], rna_common[cell_idx,j])
  
  xgb_fit <- xgboost::xgboost(data = dat_mat, 
                              label = label_vec, 
                              max.depth = 4, 
                              nround = 5, 
                              objective = "binary:logistic",
                              verbose = 0)
  
  xgb_pred <- as.numeric(xgboost:::predict.xgb.Booster(
    xgb_fit, 
    newdata = dat_mat
  ) > 0.5)
  
  confusion_mat <- table(label_vec, xgb_pred)
  
  # if(ncol(confusion_mat) == 1){
  #   val1 <- confusion_mat[1,1]/sum(confusion_mat[1,])
  #   val2 <- 0
  # } else {
  #   val1 <- confusion_mat[1,1]/sum(confusion_mat[1,])
  #   val2 <- confusion_mat[2,2]/sum(confusion_mat[2,])
  # }
  # (val1+val2)/2
  
  if(ncol(confusion_mat) == 1){
    return(0)
  } else {
    return(confusion_mat[2,2]/sum(confusion_mat[2,]))
  }
})
quantile(classification_vec)
length(which(classification_vec>0.075))

gene_idx <- order(classification_vec, decreasing = T)[1:25]

set.seed(10)
rng_idx <- sample(1:length(c(cell_idx, cell_idx2)))

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_COCL2_RNA-chromAct_leafplot_withweek5_genes_random-forest.png"),
    height = 3000, width = 3000, units = "px", res = 300)
par(mar = c(4,4,4,0.5), mfrow = c(5,5))
for(j in gene_idx){
  plot(x = atac_pred[c(cell_idx, cell_idx2)[rng_idx],j],
       y = rna_common[c(cell_idx, cell_idx2)[rng_idx],j], 
       pch = 16, col = c(color_vec, color_vec2)[rng_idx],
       xlab = "Chromatin act. value", ylab = "RNA value", 
       main = paste0(colnames(rna_common)[j], "\nCorr: ", round(cor_vec[j],2)),
       cex.main = 0.7)
}
graphics.off()
