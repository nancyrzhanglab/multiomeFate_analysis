rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

gene_vec <- intersect(gene_vec, rownames(all_data[["Saver"]]@data))
gene_vec <- sort(gene_vec)
length(gene_vec)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
treatment_vec <- c("CIS", "COCL2", "DABTRAM")

############################

treatment <- "CIS"
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 5),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 25))]
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

keep_vec <- rep(NA, ncol(all_data))
keep_vec[tier1_idx] <- paste0("3high_winner_", treatment)
keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
keep_vec[tier3_idx] <- paste0("1loser_", treatment)
table(keep_vec)
all_data2 <- all_data
all_data2$keep <- keep_vec
all_data2 <- subset(all_data2, keep %in% c(paste0("3high_winner_", treatment),
                                           paste0("2mid_winner_", treatment),
                                           paste0("1loser_", treatment)))

rna_mat <- all_data2[["Saver"]]@data[gene_vec,]
rna_mat <- t(rna_mat)
tier_vec <- all_data2$keep
dim(rna_mat); length(tier_vec); table(tier_vec)

rm(list=c("all_data", "all_data2"))
gc()

#################################################

load("../../../../out/kevin/Writeup6g/custom_peakAggregation_tmp.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

len <- sapply(result_list, length)
if(any(len == 0)) result_list <- result_list[-which(len == 0)]
gene_vec <- names(result_list)
all(colnames(rna_mat) %in% gene_vec)
if(any(!colnames(rna_mat) %in% gene_vec)){
  gene_vec <- sort(unique(intersect(names(result_list), colnames(rna_mat))))
  rna_mat <- rna_mat[,gene_vec]
  stopifnot(all(colnames(rna_mat) %in% names(result_list)))
  gene_vec <- colnames(rna_mat)
}

# result_list <- result_list[which(gene_vec %in% rownames(rna_mat))]
# result_list <- result_list[rownames(rna_mat)]

chract_mat <- sapply(result_list, function(mat){
  sparseMatrixStats::rowSums2(mat)
})
rownames(chract_mat) <- rownames(result_list[[1]])
# quantile(chract_mat); length(which(chract_mat == 0))/prod(dim(chract_mat))
# quantile(chract_mat[chract_mat != 0])

lib_mat <- Matrix::rowSums(chract_mat)
quantile(lib_mat)
lib_mat <- pmax(lib_mat, 1)

######################

# for diag(vec) %*% mat
.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

# for mat %*% diag(vec)
# see https://stackoverflow.com/questions/17080099/fastest-way-to-multiply-matrix-columns-with-vector-elements-in-r
.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}

######################

chract_mat_normalized <- .mult_vec_mat(1/lib_mat, chract_mat)
chract_mat_normalized2 <- log1p(chract_mat_normalized*1e6)
chract_mat_normalized3 <- scale(chract_mat_normalized2)

set.seed(10)
dimred_res <- irlba::irlba(chract_mat_normalized3, nv = 50)
eigenbasis <- dimred_res$u[,2:50]
chract_mat_normalized4 <- tcrossprod(eigenbasis %*% diag(dimred_res$d[2:50]), dimred_res$v[,2:50])
rownames(chract_mat_normalized4) <- rownames(chract_mat)
colnames(chract_mat_normalized4) <- colnames(chract_mat)

chr_peak_list <- result_list[colnames(rna_mat)]

for(i in 1:length(chr_peak_list)){
  print(i)
  gene <- names(chr_peak_list)[i]
  mat <- chr_peak_list[[i]]
  
  # apply library size
  mat <- .mult_vec_mat(1/lib_mat, mat)
  
  # rescale based on the log-version
  mat2 <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
  base <- chract_mat_normalized2[,gene]
  denom <- chract_mat_normalized[,gene]
  denom[base == 0] <- 1
  stopifnot(min(denom) > 0)
  for(j in 1:ncol(mat2)){
    mat2[,j] <- base * mat[,j]/denom
  }
  
  mean_val <- mean(chract_mat_normalized2[,gene])
  sd_val <- sd(chract_mat_normalized2[,gene])
  mat3 <- (mat2 - mean_val)/sd_val
  
  mat4 <- eigenbasis %*% crossprod(eigenbasis, mat3)
  rownames(mat4) <- rownames(chr_peak_list[[i]])
  colnames(mat4) <- colnames(chr_peak_list[[i]])
  
  chr_peak_list[[i]] <- mat4
}

# now trim to only the day10s
for(i in 1:length(chr_peak_list)){
  chr_peak_list[[i]] <- chr_peak_list[[i]][rownames(rna_mat),]
}

ls_vec <- ls()
ls_vec <- ls_vec[-which(ls_vec %in% c("rna_mat", "chr_peak_list", "tier_vec", "date_of_run", "session_info"))]
rm(list = ls_vec); gc()

#######################

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))
gene_vec <- sort(intersect(gene_vec, colnames(rna_mat)))

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

################

treatment <- "CIS"

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
  
  png(paste0("../../../../out/figures/Writeup6g/Writeup6g_day10_spca_leafplot_", treatment, "_", kk, ".png"),
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



