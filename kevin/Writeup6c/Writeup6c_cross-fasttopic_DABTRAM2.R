rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}

###############################

gene_vec <- intersect(rownames(all_data[["fasttopic_DABTRAM"]]@feature.loadings), 
                      rownames(all_data[["geneActivity"]]@scale.data))
cell_idx <- which(all_data$dataset %in% c("day0", "day10_DABTRAM", "week5_DABTRAM"))

mat_2 <- all_data[["geneActivity"]]@counts[gene_vec,cell_idx]
mat_2 <- .mult_mat_vec(mat_2, 1/Matrix::colSums(mat_2))
mat_2 <- t(as.matrix(mat_2))

# construct proxy scores
feature_rna <- all_data[["fasttopic_DABTRAM"]]@feature.loadings[gene_vec,]
proxy_atac <- mat_2 %*% feature_rna
proxy_atac <- proxy_atac/median(proxy_atac)
# for(i in 1:nrow(proxy_atac)){
#   proxy_atac[i,] <- proxy_atac[i,]/sum(proxy_atac[i,])
# }
rownames(proxy_atac) <- colnames(all_data)[cell_idx]

#############################################################

# p-value
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

atac_p_value_vec <- sapply(1:ncol(proxy_atac), function(j){
  res <- stats::wilcox.test(x = proxy_atac[cell_winning,j], 
                            y = proxy_atac[cell_losing,j])
  res$p.value
})
rna_p_value_vec <- sapply(1:ncol(rna_loading), function(j){
  res <- stats::wilcox.test(x = rna_loading[cell_winning,j], 
                            y = rna_loading[cell_losing,j])
  res$p.value
})

cell_day0_idx <- intersect(which(all_data$dataset == "day0"),
                           which(all_data$assigned_posterior > 0.5))
cell_week5_idx <- intersect(which(all_data$dataset == "week5_DABTRAM"),
                            which(all_data$assigned_posterior > 0.5))
cell_day0 <-  colnames(all_data)[cell_day0_idx]
cell_week5 <-  colnames(all_data)[cell_week5_idx]

# Chromatin Activity plots
p <- vector("list",30)
for(j in 1:ncol(proxy_atac)){
  
  Condition <- c(rep("D0", length(cell_day0)),
                 rep("D10-Loss", length(cell_losing)),
                 rep("D10-Win", length(cell_winning)),
                 rep("W5", length(cell_week5)))
  Embedding <- c(proxy_atac[c(cell_day0,
                              cell_losing,
                              cell_winning,
                              cell_week5),j])
  tab <- data.frame(Condition = Condition, 
                    Embedding = Embedding)
  
  p[[j]] <- ggplot2::ggplot(tab, ggplot2::aes(x=Condition, y=Embedding, fill=Condition))
  p[[j]] <- p[[j]] + ggplot2::geom_violin(trim=FALSE, show.legend=FALSE, scale = "width")
  p[[j]] <- p[[j]] + ggplot2::geom_boxplot(width=0.1, show.legend=FALSE)
  p[[j]] <- p[[j]] + ggplot2::ggtitle(paste("Topic",j, ", -Log10(p)=",round(-log10(atac_p_value_vec[j]),2)))
}

pdf(file = "../../../../out/figures/Writeup6c/Writeup6c_fasttopic_DABTRAM_chromatinActivity2.pdf", 
    width = 18, height = 14)
cowplot::plot_grid(plotlist = p, ncol = 6, nrow = 5)
dev.off()

###################################

color_palette <- scales::hue_pal()(4)
color_vec <- c(rep(color_palette[1], length(cell_day0)),
               rep(color_palette[2], length(cell_losing)),
               rep(color_palette[3], length(cell_winning)),
               rep(color_palette[4], length(cell_week5)))
important_topics <- c(1,7,11,12,13,17,25,27,30)
rna_loading <- all_data[["fasttopic_DABTRAM"]]@cell.embeddings

set.seed(10)
rng_idx <- sample(1:length(color_vec))

for(k in important_topics){
  print(k)
  png(file = paste0("../../../../out/figures/Writeup6c/Writeup6c_fasttopic_DABTRAM_RNA-chromatinActivity_topic", k, ".png"), 
      width = 1500, height = 1500, res = 300, units = "px")
  vec1 <- c(proxy_atac[c(cell_day0,
                         cell_losing,
                         cell_winning,
                         cell_week5),k])
  vec2 <- c(rna_loading[c(cell_day0,
                          cell_losing,
                          cell_winning,
                          cell_week5),k])
  cor_val <- stats::cor(vec1, vec2, method = "kendall")
  
  xlim <- c(0, quantile(vec1, probs = 0.95))
  ylim <- c(0, quantile(vec2, probs = 0.95))
  
  plot(x = vec1[rng_idx], y = vec2[rng_idx], 
       xlab = "Chromatin Activity", 
       ylab = "RNA",
       main = paste0("Topic ", k, ", Corr: ", round(cor_val, 2)),
       pch = 16, col = color_vec[rng_idx],
       xlim = xlim, ylim = ylim)
  graphics.off()
}

kendall_vec <- sapply(1:30, function(k){
  print(k)
  vec1 <- c(proxy_atac[c(cell_day0,
                         cell_losing,
                         cell_winning,
                         cell_week5),k])
  vec2 <- c(rna_loading[c(cell_day0,
                          cell_losing,
                          cell_winning,
                          cell_week5),k])
  stats::cor(vec1, vec2, method = "kendall")
})

png(file = paste0("../../../../out/figures/Writeup6c/Writeup6c_fasttopic_DABTRAM_RNA-chromatinActivity_kendall.png"), 
    width = 2500, height = 1000, res = 300, units = "px")
barplot(kendall_vec, main="Correlation (Kendall) between\nRNA and Chromatin Activity loading",
        xlab="Topic #")
graphics.off()

####################################################3

# make dedicated violin plots

Condition <- c(rep("D0", length(cell_day0)),
               rep("D10-Loss", length(cell_losing)),
               rep("D10-Win", length(cell_winning)),
               rep("W5", length(cell_week5)))

for(topic_num in important_topics){
  p <- vector("list",2)
  for(j in 1:2){
    if(j == 1){
      Embedding <- c(rna_loading[c(cell_day0,
                                   cell_losing,
                                   cell_winning,
                                   cell_week5),topic_num])
    } else {
      Embedding <- c(proxy_atac[c(cell_day0,
                                  cell_losing,
                                  cell_winning,
                                  cell_week5),topic_num])
    }
    
    tab <- data.frame(Condition = Condition, 
                      Embedding = Embedding)
    
    p[[j]] <- ggplot2::ggplot(tab, ggplot2::aes(x=Condition, y=Embedding, fill=Condition))
    p[[j]] <- p[[j]] + ggplot2::geom_violin(trim=FALSE, show.legend=FALSE, scale = "width")
    p[[j]] <- p[[j]] + ggplot2::geom_boxplot(width=0.1, show.legend=FALSE)
    if(j == 1){
      topic_name <- "RNA, Topic "; pval <- rna_p_value_vec[topic_num]
    } else {
      topic_name <- "ChrAct, Topic "; pval <- atac_p_value_vec[topic_num]
    }
    p[[j]] <- p[[j]] + ggplot2::ggtitle(paste(topic_name,topic_num, ", -Log10(p)=",round(-log10(pval),2)))
  }
  
  pdf(file = paste0("../../../../out/figures/Writeup6c/Writeup6c_fasttopic_DABTRAM_topic_", topic_num, ".pdf"), 
      width = 10, height = 5)
  cowplot::plot_grid(plotlist = p, ncol = 2, nrow = 1)
  dev.off()
}

##################################

feature_rna2 <- feature_rna
for(k in 1:ncol(feature_rna2)){
  feature_rna2[,k] <- feature_rna2[,k]/sum(feature_rna2[,k])
}

# now plot the weights of the genes
for(k in important_topics){
  png(file = paste0("../../../../out/figures/Writeup6c/Writeup6c_fasttopic_DABTRAM_gene-weights_topic", k, ".png"), 
      width = 1200, height = 1000, res = 300, units = "px")
  tmp <- hist(feature_rna2[,k], breaks = 21)
  tmp$counts <- log10(tmp$counts + 1)
  plot(tmp,
       main = paste0("Gene weights for Topic ", k),
       xlab = "Gene weights (Sums to 1)",
       ylab = "Log10 of frequency",
       col = "gray")
  mean_val <- mean(feature_rna2[,k])
  median_val <- median(feature_rna2[,k])
  
  lines(rep(mean_val, 2), c(-1e7,1e7), col = 2, lwd = 2)
  lines(rep(median_val, 2), c(-1e7,1e7), col = 3, lwd = 2, lty = 2)
  lines(c(-100,100), rep(log10(1+1), 2), lty = 3)
  graphics.off()
}

######################

k <- 1
round(sort(feature_rna2[,k], decreasing = T)[1:10],3)

k <- 7
round(sort(feature_rna2[,k], decreasing = T)[1:20],3)
zz <- sort(feature_rna2[,k], decreasing = T); zz <- zz[zz>=0.0095]
names(zz); length(zz)

k <- 11
round(sort(feature_rna2[,k], decreasing = T)[1:20],3)
zz <- sort(feature_rna2[,k], decreasing = T); zz <- zz[zz>=0.0095]
names(zz); length(zz)

k <- 25
round(sort(feature_rna2[,k], decreasing = T)[1:20],3)
zz <- sort(feature_rna2[,k], decreasing = T); zz <- zz[zz>=0.0095]
names(zz); length(zz)

k <- 27
round(sort(feature_rna2[,k], decreasing = T)[1:20],3)
zz <- sort(feature_rna2[,k], decreasing = T); zz <- zz[zz>=0.0095]
names(zz); length(zz)

k <- 30
round(sort(feature_rna2[,k], decreasing = T)[1:20],3)
zz <- sort(feature_rna2[,k], decreasing = T); zz <- zz[zz>=0.0095]
names(zz); length(zz)

important_topics2 <- c(1,7,11,25,27,30)
round(cor(rna_loading[cell_idx,])[important_topics2,important_topics2],1)
round(cor(feature_rna)[important_topics2,important_topics2],1)


