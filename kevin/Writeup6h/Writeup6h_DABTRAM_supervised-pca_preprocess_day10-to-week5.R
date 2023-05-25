rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)
load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM_supervised-pca_preprocess.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

###########

chr_peak_list <- preprocess_res$chr_peak_list
rna_mat <- preprocess_res$rna_mat
gene_vec <- colnames(rna_mat)

treatment <- "DABTRAM"
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]

tier1_idx <- intersect(
  intersect(which(metadata_mat$assigned_lineage %in% tier1_lineages),
            which(metadata_mat$assigned_posterior >= 0.5)),
  which(metadata_mat$dataset == paste0("day10_", treatment))
)
tier2_idx <- intersect(
  intersect(which(metadata_mat$assigned_lineage %in% tier2_lineages),
            which(metadata_mat$assigned_posterior >= 0.5)),
  which(metadata_mat$dataset == paste0("day10_", treatment))
)
tier3_idx <- intersect(
  intersect(which(metadata_mat$assigned_lineage %in% tier3_lineages),
            which(metadata_mat$assigned_posterior >= 0.5)),
  which(metadata_mat$dataset == paste0("day10_", treatment))
)

tier_vec <- rep(NA, nrow(metadata_mat))
names(tier_vec) <- rownames(metadata_mat)
tier_vec[tier1_idx] <- paste0("3high_winner_", treatment)
tier_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
tier_vec[tier3_idx] <- paste0("1loser_", treatment)
tier_vec <- tier_vec[!is.na(tier_vec)]

rna_mat <- rna_mat[names(tier_vec),,drop = F]
chr_peak_list <- lapply(chr_peak_list, function(mat){
  mat[names(tier_vec),,drop = F]
})

##################

cell_names <- names(tier_vec)[which(!is.na(tier_vec))]

tier_vec <- tier_vec[rownames(rna_mat)]
y <- multiomeFate:::form_onehot_classification_mat(tier_vec)

spca_res_list <- vector("list", length = length(gene_vec))
names(spca_res_list) <- gene_vec
for(i in 1:length(gene_vec)){
  gene <- gene_vec[i]
  print(paste0(gene, ": ", i, " out of ", length(gene_vec)))
  if(stats::sd(rna_mat[,gene]) <= 1e-6) next()
  
  chr_peak_mat <- chr_peak_list[[gene]]
  if(!is.matrix(chr_peak_mat)) {
    chr_peak_mat <- matrix(chr_peak_mat, nrow = length(chr_peak_mat), ncol = 1)
    colnames(chr_peak_mat) <- paste0(gene, ":ATAC")
  }
  tmp <- cbind(rna_mat[,gene], chr_peak_mat)
  colnames(tmp)[1] <- paste0(gene, ":RNA")
  sd_vec <- apply(tmp, 2, stats::sd)
  if(any(sd_vec <= 1e-6)){
    tmp <- tmp[,sd_vec >= 1e-6,drop=F]
  }
  tmp <- scale(tmp)
  spca_res_list[[gene]] <- multiomeFate:::supervised_pca(x = tmp, y = y,
                                                         k = 2,
                                                         scale_x = T,
                                                         orthogonalize = T)
}

save(spca_res_list, 
     tab_mat, treatment, tier_vec,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM_supervised-pca_day10-to-week5.RData")

##########

source("../Writeup6g/ordinal_functions.R")

y_vec <- rep(NA, length(tier_vec))
idx_list <- vector("list", length = 3)
idx_list[[1]] <- which(tier_vec == paste0("1loser_", treatment))
idx_list[[2]] <- which(tier_vec == paste0("2mid_winner_", treatment))
idx_list[[3]] <- which(tier_vec == paste0("3high_winner_", treatment))
for(i in 1:3){
  y_vec[idx_list[[i]]] <- i
}
y_vec <- as.factor(y_vec)

cv_score_vec <- rep(NA, length(gene_vec))
names(cv_score_vec) <- gene_vec
for(i in 1:length(gene_vec)){
  gene <- gene_vec[i]
  print(paste0(gene, ": ", i, " out of ", length(gene_vec)))
  set.seed(10)
  
  x_mat <- spca_res_list[[gene]]$dimred
  
  if(length(x_mat) == 0) next()
  cv_score_vec[gene] <- .five_fold_cv(x_mat, y_vec)
}

save(spca_res_list, cv_score_vec,
     tab_mat, treatment, tier_vec,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM_supervised-pca_day10-to-week5.RData")

