rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)
library(ordinal)

load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak.RData")
treatment <- "DABTRAM"

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

gene_vec <- intersect(gene_vec, colnames(rna_mat))
gene_vec <- sort(gene_vec)
length(gene_vec)

n <- nrow(rna_mat)
y <- multiomeFate:::form_onehot_classification_mat(tier_vec)

spca_res_list <- vector("list", length = length(gene_vec))
names(spca_res_list) <- gene_vec
for(i in 1:length(gene_vec)){
  gene <- gene_vec[i]
  print(paste0(gene, ": ", i, " out of ", length(gene_vec)))
  tmp <- cbind(rna_mat[,gene], chr_peak_list[[gene]])
  colnames(tmp)[1] <- paste0(gene, ":RNA")
  tmp <- scale(tmp)
  spca_res_list[[gene]] <- multiomeFate:::supervised_pca(x = tmp, y = y)
}

# head(spca_res_list[["FN1"]]$U)
# head(spca_res_list[["FN1"]]$dimred)
# 
# percent_mat <- sapply(spca_res_list, function(zz){
#   tmp <- zz$U
#   (tmp[1,]^2)*100
# })
# round(percent_mat)

##############################

y_vec <- rep(NA, length(tier_vec))
idx_list <- vector("list", length = 3)
idx_list[[1]] <- which(tier_vec == paste0("1loser_", treatment))
idx_list[[2]] <- which(tier_vec == paste0("2mid_winner_", treatment))
idx_list[[3]] <- which(tier_vec == paste0("3high_winner_", treatment))
for(i in 1:3){
  y_vec[idx_list[[i]]] <- i
}
y_vec <- as.factor(y_vec)

cv_score_vec <- rep(NA, length(spca_res_list))
names(cv_score_vec) <- names(spca_res_list)
for(gene in names(spca_res_list)){
  print(gene)
  set.seed(10)
  
  x_mat <- spca_res_list[[gene]]$dimred
  cv_score_vec[gene] <- .five_fold_cv(x_mat, y_vec)
}
round(quantile(100*cv_score_vec))




# y_vec <- rep(NA, length(tier_vec))
# weight_vec <- rep(NA, length(tier_vec))
# idx_list <- vector("list", length = 3)
# idx_list[[1]] <- which(tier_vec == paste0("1loser_", treatment))
# idx_list[[2]] <- which(tier_vec == paste0("2mid_winner_", treatment))
# idx_list[[3]] <- which(tier_vec == paste0("3high_winner_", treatment))
# len_vec <- sapply(idx_list, length)
# 
# for(i in 1:3){
#   y_vec[idx_list[[i]]] <- i
#   weight_vec[idx_list[[i]]] <- 1/length(idx_list[[i]])
# }
# 
# gene <- "FN1"
# tmp_df <- cbind(y_vec, spca_res_list[[gene]]$dimred)
# colnames(tmp_df)[1] <- "y"
# tmp_df <- as.data.frame(tmp_df)
# tmp_df[,"y"] <- as.factor(tmp_df[,"y"])
# 
# oridinal_res <- ordinal::clm(y ~ ., data = tmp_df, weights = weight_vec)
# # summary(oridinal_res)
# pred <- stats::predict(oridinal_res, newdata = subset(tmp_df, select = -y))$fit
# pred <- sapply(1:nrow(pred), function(i){which.max(pred[i,])})
# # table(pred)





