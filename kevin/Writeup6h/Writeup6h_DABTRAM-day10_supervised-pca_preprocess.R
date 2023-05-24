rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6h/peak-gene-matching.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
metadata <- all_data@meta.data

treatment <- "DABTRAM"
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)

tier1_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier1_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier2_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier2_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier3_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier3_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)

tier_vec <- rep(NA, ncol(all_data))
names(tier_vec) <- colnames(all_data)
tier_vec[tier1_idx] <- paste0("3high_winner_", treatment)
tier_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
tier_vec[tier3_idx] <- paste0("1loser_", treatment)

all_data2 <- all_data
all_data2$tier <- tier_vec
all_data2 <- subset(all_data2, tier %in% c(paste0("3high_winner_", treatment),
                                           paste0("2mid_winner_", treatment),
                                           paste0("1loser_", treatment)))

preprocess_res <- extract_relevant_peaks(peak_mapping_list = matching_list,
                                            seurat_obj = all_data2,
                                            slot_atac = "ATAC",
                                            slot_rna = "Saver",
                                            verbose = 3)
preprocess_res <- preprocess_chromatin_peak(chr_peak_list = preprocess_res$chr_peak_list,
                                            rna_mat = preprocess_res$rna_mat,
                                            verbose = 3)

save(preprocess_res, tab_mat, treatment, tier_vec,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_supervised-pca_preprocess.RData")

#######

rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)
load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_supervised-pca_preprocess.RData")

chr_peak_list <- preprocess_res$chr_peak_list
rna_mat <- preprocess_res$rna_mat

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))
gene_vec <- intersect(gene_vec, colnames(rna_mat))


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

n <- nrow(rna_mat)
col_vec <- rep(NA, n)
color_palette <- grDevices::colorRampPalette(c("lightgray", "dodgerblue4"))(3)
idx1 <- which(tier_vec ==  paste0("1loser_", treatment))
idx2 <- which(tier_vec ==  paste0("2mid_winner_", treatment))
idx3 <- which(tier_vec ==  paste0("3high_winner_", treatment))
idx_list <- list(idx1,idx2,idx3)
col_vec[idx1] <- color_palette[1]
col_vec[idx2] <- color_palette[2]
col_vec[idx3] <- color_palette[3]
order_idx <- c(idx1,idx2,idx3)

for(kk in 1:2){
  gene_vec_tmp <- gene_vec[((kk-1)*25+1):min(kk*25, length(gene_vec))]
  
  png(paste0("../../../../out/figures/Writeup6h/Writeup6h_day10_spca_leafplot_", treatment, "_", kk, ".png"),
      height = 3000, width = 2500, res = 300, units = "px")
  par(mfrow = c(5,5), mar = c(4,4,3,0.5))
  for(gene in gene_vec_tmp){
    tmp <-  Re(spca_res_list[[gene]]$dimred)
    vec1 <- tmp[,1]
    vec2 <- tmp[,2]
    xlim <- quantile(vec1, probs = c(0.01, 0.99))
    ylim <- quantile(vec2, probs = c(0.01, 0.99))
    
    plot(x = vec1[order_idx], y = vec2[order_idx],
         xlim = xlim, ylim = ylim, 
         pch = 16, col = col_vec[order_idx],
         xlab = "Supervised PC1", ylab = "Supervised PC2", 
         main = gene)
    
    mean_mat <- sapply(1:3, function(i){
      matrixStats::colMedians(tmp[idx_list[[i]],])
    })
    for(i in 1:3){
      points(mean_mat[1,i], mean_mat[2,i], 
             col = "black", pch = 16, cex = 4)
      points(mean_mat[1,i], mean_mat[2,i], 
             col = "white", pch = 16, cex = 3.4)
    } 
    for(i in 1:3){
      points(mean_mat[1,i], mean_mat[2,i], 
             col = color_palette[i], pch = 16, cex = 3)
    }
    
  }
  graphics.off()
}

