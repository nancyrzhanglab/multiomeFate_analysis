rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)
load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM_supervised-pca_day0-to-day10.RData")

gene_vec <- names(spca_res_list)
spca_rna_importance <- sapply(spca_res_list, function(lis){
  if(length(lis) == 0) return(NULL)
  
  var_mat <- stats::var(lis$dimred)
  var_vec <- diag(var_mat)
  rna_percentage <- sapply(1:ncol(lis$U), function(j){
    lis$U[1,j]^2
  })
  
  (rna_percentage %*% var_vec)/sum(var_vec)
})
idx <- unique(c(which(sapply(spca_rna_importance, is.null)), 
                which(sapply(cv_score_vec, is.null))))
if(length(idx) != 0){
  spca_rna_importance <- unlist(spca_rna_importance[-idx])
  cv_score_vec <- cv_score_vec[-idx]
  gene_vec <- names(spca_rna_importance)
  stopifnot(all(names(spca_rna_importance) == names(cv_score_vec)))
}

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec2 <- sort(unique(c(unlist(keygenes), keygenes_csc)))

df <- data.frame(cv_score_vec = cv_score_vec,
                 gene = gene_vec,
                 labeling = gene_vec %in% gene_vec2,
                 spca_rna_importance = spca_rna_importance)
# put all the labeling == TRUE on bottom
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = spca_rna_importance, y = cv_score_vec))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = gene, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM_supervised-pca_day0-to-day10_all-genes.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

#######################################

gene_vec2 <- names(spca_rna_importance)[spca_rna_importance >= 0.5]

df <- data.frame(cv_score_vec = cv_score_vec,
                 gene = gene_vec,
                 labeling = gene_vec %in% gene_vec2,
                 spca_rna_importance = spca_rna_importance)
# put all the labeling == TRUE on bottom
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = spca_rna_importance, y = cv_score_vec))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = gene, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM_supervised-pca_day0-to-day10_high-rna.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

##############################

tier_vec_save <- tier_vec
load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM_supervised-pca_preprocess.RData")

chr_peak_list <- preprocess_res$chr_peak_list
rna_mat <- preprocess_res$rna_mat
tier_vec <- tier_vec_save

rna_mat <- rna_mat[names(tier_vec),,drop = F]
chr_peak_list <- lapply(chr_peak_list, function(mat){
  mat[names(tier_vec),,drop = F]
})

gene_vec2 <- names(cv_score_vec)[order(cv_score_vec, decreasing = T)[1:50]]

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
  gene_vec_tmp <- gene_vec2[((kk-1)*25+1):min(kk*25, length(gene_vec2))]
  
  png(paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM_supervised-pca_day0-to-day10_leafplot_", treatment, "_largest-cv_", kk, ".png"),
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
         main = paste0(gene, " (CV: ", round(cv_score_vec[gene],2), ")"),
         cex.main = 0.9)
    
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
