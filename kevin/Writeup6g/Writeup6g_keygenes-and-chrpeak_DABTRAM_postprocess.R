rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)
library(ordinal)

load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM.RData")
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

########

# null_dist_list <- vector("list", length = length(spca_res_list))
# names(null_dist_list) <- names(spca_res_list)
# for(gene in names(spca_res_list)){
#   print(gene)
#   x_mat <- spca_res_list[[gene]]$dimred
#   
#   set.seed(10)
#   null_dist_list[[gene]] <- .permutation_null_score(x_mat, y_vec, trials = 100)
# }
# 
# p_val <- sapply(1:length(null_dist_list), function(i){
#   length(which(cv_score_vec[i] <= null_dist_list[[i]]))/length(null_dist_list[[i]])
# })

##########

percent_rna <- sapply(gene_vec, function(gene){
  100*spca_res_list[[gene]]$U[1,1]^2
})

df <- data.frame(cv_score_vec = cv_score_vec,
                 gene = gene_vec,
                 percent_rna = percent_rna)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = percent_rna, y = cv_score_vec))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggrepel::geom_text_repel(data = df,
                                    ggplot2::aes(label = gene),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_rnapercentage_cvscore.png"),
                p1, device = "png", width = 10, height = 10, units = "in")
                                                                                            