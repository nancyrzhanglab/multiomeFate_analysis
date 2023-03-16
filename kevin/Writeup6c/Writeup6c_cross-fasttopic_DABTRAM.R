rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6b/Writeup6b_chromatinAct_fasttopics_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

rna_loading <- all_data[["fasttopic_DABTRAM"]]@feature.loadings
atac_loading <- topic_res$F

gene_vec <- intersect(rownames(rna_loading), rownames(atac_loading))
cor_mat <- cor(sqrt(rna_loading[gene_vec,]), sqrt(atac_loading[gene_vec,]))
# not that impressive...

#################################

# let's do a janky transformation
# initialize
cell_idx <- which(all_data$dataset %in% c("day0", "day10_DABTRAM", "week5_DABTRAM"))
mat_2 <- all_data[["geneActivity"]]@counts[gene_vec,cell_idx]
mat_2 <- t(as.matrix(mat_2))
feature_2 <- all_data[["fasttopic_DABTRAM"]]@feature.loadings[gene_vec,]
feature_sums <- colSums(feature_2)
feature_2 <- feature_2 %*% diag(1/feature_sums)
cellsum <- rowSums(mat_2)
loading_old <- matrix(cellsum, nrow = nrow(mat_2), ncol = ncol(feature_2))

###################################
# loading_old <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[cell_idx,]
# rna_colsum <- Matrix::colSums(all_data[["RNA"]]@counts[gene_vec,cell_idx])
# loading_old <- diag(rna_colsum) %*% loading_old
#
# # do an adjustment so loading and features are on the same scale
# loading_l2 <- apply(loading_old, 2, function(x){sqrt(sum(x^2))})
# feature_l2 <- apply(feature_2, 2, function(x){sqrt(sum(x^2))})
# feature_adjustment <- sqrt(loading_l2)/sqrt(feature_l2)
# loading_adjustment <- sqrt(feature_l2)/sqrt(loading_l2)
# feature_2 <- feature_2 %*% diag(feature_adjustment)
# loading_old <- loading_old %*% diag(loading_adjustment)
# 
# # add a bit of randomness so correlation isn't too high
# max_vec <- as.numeric(loading_old)/2
# set.seed(10)
# rng_vec <- stats::runif(length(max_vec), min = 0, max = max_vec)
# loading_old <- loading_old + matrix(rng_vec, nrow = nrow(loading_old), ncol = ncol(loading_old))
###################################

# from https://github.com/stephenslab/fastTopics/blob/master/R/misc.R
scale.cols <- function (A, b)
  t(t(A) * b)

# from https://github.com/stephenslab/fastTopics/blob/master/R/betanmf.R
betanmf_update_loadings <- function (X, A, B)
  scale.cols(A * tcrossprod(X / (A %*% B),B),1/rowSums(B))

loading_new <- betanmf_update_loadings(X = mat_2, 
                                       A = loading_old,
                                       B = t(feature_2))
sum(abs(mat_2 - tcrossprod(loading_new,feature_2)))/prod(dim(mat_2))

# reparameterize the cell loadings)
feature_sums <- colSums(feature_2)
size_vec <- as.numeric(loading_new %*% feature_sums)
loading_new2 <- diag(1/size_vec) %*% loading_new %*% diag(feature_sums)
quantile(rowSums(loading_new2))
rownames(loading_new2) <- rownames(loading_new)
colnames(loading_new2) <- colnames(loading_new)

##########

loading_rna <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings[cell_idx,]
colnames(loading_new2) <- paste0("ChrAct_", colnames(loading_new2))
cor_mat <- cor(loading_new2, loading_rna)
important_topics <- c(1,5,11,12,13,17,27,30)
round(cor_mat[,important_topics], 2)
zz <- cbind(round(apply(cor_mat, 2, max),2), 
            apply(cor_mat, 2, which.max),
            round(apply(cor_mat, 2, function(x){sort(x, decreasing = T)[2]}),2),
            round(diag(cor_mat),2))
colnames(zz) <- c("Max", "Maximizing idx", "Second max", "Diagonal")
zz
zz[important_topics,]

cor_mat2 <- cor(feature_2)
round(cor_mat2,1)
cor_mat2[lower.tri(cor_mat2,diag = T)] <- 0
which(cor_mat2 >= 0.6, arr.ind = T)

###################################

# compute the p-values
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_DABTRAM"] >= 50)]
cell_day10_idx <- intersect(which(all_data$dataset == "day10_DABTRAM"),
                            which(all_data$assigned_posterior > 0.5))
cell_winning_idx <- intersect(cell_day10_idx, 
                              which(all_data$assigned_lineage %in% lineage_vec))
cell_losing_idx <- intersect(cell_day10_idx, 
                             which(!all_data$assigned_lineage %in% lineage_vec))
cell_winning <- colnames(all_data)[cell_winning_idx]
cell_losing <- colnames(all_data)[cell_losing_idx]

atac_p_value_vec <- sapply(1:ncol(loading_new2), function(j){
  print(j)
  res <- stats::wilcox.test(x = loading_new2[cell_winning,j], 
                            y = loading_new2[cell_losing,j])
  res$p.value
})
rna_p_value_vec <- sapply(1:ncol(loading_rna), function(j){
  print(j)
  res <- stats::wilcox.test(x = loading_rna[cell_winning,j], 
                            y = loading_rna[cell_losing,j])
  res$p.value
})

# make some plots
cell_day0_idx <- intersect(which(all_data$dataset == "day0"),
                           which(all_data$assigned_posterior > 0.5))
cell_week5_idx <- intersect(which(all_data$dataset == "week5_DABTRAM"),
                            which(all_data$assigned_posterior > 0.5))
cell_day0 <-  colnames(all_data)[cell_day0_idx]
cell_week5 <-  colnames(all_data)[cell_week5_idx]

# RNA plots
p <- vector("list",30)
for(j in 1:ncol(loading_rna)){
  
  Condition <- c(rep("D0", length(cell_day0)),
                 rep("D10-Loss", length(cell_losing)),
                 rep("D10-Win", length(cell_winning)),
                 rep("W5", length(cell_week5)))
  Embedding <- c(loading_rna[c(cell_day0,
                               cell_losing,
                               cell_winning,
                               cell_week5),j])
  tab <- data.frame(Condition = Condition, 
                    Embedding = Embedding)
  
  p[[j]] <- ggplot2::ggplot(tab, ggplot2::aes(x=Condition, y=Embedding, fill=Condition))
  p[[j]] <- p[[j]] + ggplot2::geom_violin(trim=FALSE, show.legend=FALSE, scale = "width")
  p[[j]] <- p[[j]] + ggplot2::geom_boxplot(width=0.1, show.legend=FALSE)
  p[[j]] <- p[[j]] + ggplot2::ggtitle(paste("Topic",j, ", -Log10(p)=",round(-log10(rna_p_value_vec[j]),2)))
}

pdf(file = "../../../../out/figures/Writeup6c/Writeup6c_fasttopic_DABTRAM_RNA.pdf", width = 18, height = 14)
cowplot::plot_grid(plotlist = p, ncol = 6, nrow = 5)
dev.off()

# Chromatin Activity plots
p <- vector("list",30)
for(j in 1:ncol(loading_new2)){
  
  Condition <- c(rep("D0", length(cell_day0)),
                 rep("D10-Loss", length(cell_losing)),
                 rep("D10-Win", length(cell_winning)),
                 rep("W5", length(cell_week5)))
  Embedding <- c(loading_new2[c(cell_day0,
                               cell_losing,
                               cell_winning,
                               cell_week5),j])
  tab <- data.frame(Condition = Condition, 
                    Embedding = Embedding)
  
  p[[j]] <- ggplot2::ggplot(tab, ggplot2::aes(x=Condition, y=Embedding, fill=Condition))
  p[[j]] <- p[[j]] + ggplot2::geom_violin(trim=FALSE, show.legend=FALSE)
  p[[j]] <- p[[j]] + ggplot2::geom_boxplot(width=0.1, show.legend=FALSE)
  p[[j]] <- p[[j]] + ggplot2::ggtitle(paste("Topic",j, ", -Log10(p)=",round(-log10(atac_p_value_vec[j]),2)))
}

pdf(file = "../../../../out/figures/Writeup6c/Writeup6c_fasttopic_DABTRAM_chromatinActivity.pdf", width = 18, height = 14)
cowplot::plot_grid(plotlist = p, ncol = 6, nrow = 5)
dev.off()
