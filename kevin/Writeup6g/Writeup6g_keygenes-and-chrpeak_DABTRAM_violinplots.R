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

#####################

y_vec <- rep(NA, length(tier_vec))
idx_list <- vector("list", length = 3)
idx_list[[1]] <- which(tier_vec == paste0("1loser_", treatment))
idx_list[[2]] <- which(tier_vec == paste0("2mid_winner_", treatment))
idx_list[[3]] <- which(tier_vec == paste0("3high_winner_", treatment))
for(i in 1:3){
  y_vec[idx_list[[i]]] <- i
}
y_vec <- as.factor(y_vec)

.form_df <- function(gene){
  tmp <- cbind(rna_mat[,gene], chr_peak_list[[gene]])
  colnames(tmp)[1] <- paste0(gene, ":RNA")
  tmp2 <- scale(tmp)
  
  rna_vec_org <- tmp2[,1]
  chract_vec_org <- scale(rowSums(tmp[,-1,drop=F]))
  spca_res <- multiomeFate:::supervised_pca(x = tmp2, y = y)
  
  pca_vec <- scale(spca_res$dimred[,1])
  
  tmp3 <- tcrossprod(tmp2 %*% spca_res$U, spca_res$U)
  rna_vec_denoised <- scale(tmp3[,1])
  chract_vec_denoised <- scale(rowSums(tmp3[,-1,drop=F]))
  
  df <- data.frame("RNA_original" = rna_vec_org,
                   "ChrAct_original" = chract_vec_org,
                   "RNA_SPCA" = rna_vec_denoised,
                   "ChrAct_SPCA" = chract_vec_denoised,
                   "Leading_SPCA" = pca_vec,
                   "Status" = y_vec)
  
  df
}

violin_df_list <- lapply(gene_vec, function(gene){
  print(gene)
  .form_df(gene)
})

save(violin_df_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_violinplot.RData")

#####################

rm(list=ls())
load("../../out/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_violinplot.RData")
library(ggplot2)
library(ggpubr)

df <- violin_df_list[[1]]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=RNA_original)) +
  ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
  ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::geom_boxplot(width=0.05)

p2 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=ChrAct_original)) +
  ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
  ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::geom_boxplot(width=0.05)

p3 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=RNA_SPCA)) +
  ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
  ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::geom_boxplot(width=0.05)

p4 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=ChrAct_SPCA)) +
  ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
  ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::geom_boxplot(width=0.05)

p5 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=Leading_SPCA)) +
  ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
  ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
  ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::geom_boxplot(width=0.05)

p <- ggpubr::ggarrange(p1, p3, p5, p2, p4, ncol = 3, nrow = 2) 

# ggsave("arrangedplot.png", arrange, width = 8, height = 6)

############################


pdf(paste0("../../out/figures/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_violinplot.pdf"),
    onefile = T, width = 9, height = 4.5)
for(gene in names(violin_df_list)){
  df <- violin_df_list[[gene]]
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=RNA_original)) +
    ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
    ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
    ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
    Seurat::NoLegend() + 
    ggplot2::geom_boxplot(width=0.1)
  
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=ChrAct_original)) +
    ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
    ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
    ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
    Seurat::NoLegend() + 
    ggplot2::geom_boxplot(width=0.1)
  
  p3 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=RNA_SPCA)) +
    ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
    ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
    ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
    Seurat::NoLegend() + 
    ggplot2::geom_boxplot(width=0.1)
  
  p4 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=ChrAct_SPCA)) +
    ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
    ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
    ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
    Seurat::NoLegend() + 
    ggplot2::geom_boxplot(width=0.1)
  
  p5 <- ggplot2::ggplot(df, ggplot2::aes(x=Status, y=Leading_SPCA)) +
    ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=Status)) +
    ggplot2::scale_fill_manual(values=c("lightgray", "#7190AF", "dodgerblue4")) +
    ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.3, size = 0.5) + 
    Seurat::NoLegend() + 
    ggplot2::geom_boxplot(width=0.1)
  
  p <- ggpubr::ggarrange(p1, p2, p5, p3, p4, ncol = 3, nrow = 2) 
  
  p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(gene))
  print(p)
}

dev.off() 

