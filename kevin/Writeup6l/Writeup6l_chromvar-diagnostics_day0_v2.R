rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(GGally)
library(ggplot2)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")
tab_mat_full <- table(all_data$assigned_lineage, all_data$dataset)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0_chromvar_lightweight_noATAC.RData")

tab_vec <- table(all_data$assigned_lineage)
tab_vec <- tab_vec[tab_vec > 0]
lineage_names <- names(tab_vec)
tab_mat_full <- tab_mat_full[lineage_names,]

name_vec <- names(data.use) 
motif_idx <- sort(unique(c(grep("JUN", name_vec),
                           grep("FOS", name_vec),
                           grep("NFE2", name_vec),
                           grep("TEAD", name_vec))))

# compute which cells belong to which lineage
lineage_idx_list <- lapply(lineage_names, function(lineage){
  intersect(which(all_data$assigned_lineage == lineage),
            which(all_data$dataset == "day0"))
})

# for each lineage, compute the mean scores
mean_mat <- sapply(motif_idx, function(i){
  sapply(lineage_idx_list, function(idx_vec){
    mean(all_data[["chromvar"]]@data[i,idx_vec])
  })
})
rownames(mean_mat) <- lineage_names
colnames(mean_mat) <- name_vec[motif_idx]

max_mat <- sapply(motif_idx, function(i){
  sapply(lineage_idx_list, function(idx_vec){
    max(all_data[["chromvar"]]@data[i,idx_vec])
  })
})
rownames(max_mat) <- lineage_names
colnames(max_mat) <- name_vec[motif_idx]

# make a ton of plots now
n <- nrow(tab_mat_full)
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v2-DABTRAM_ap1motif-day0size-day10size.pdf", 
    onefile = T, width = 8, height = 8)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat_full[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_dabtram = log10(tab_mat_full[,"day10_DABTRAM"] + runif(n, min = 1, max = 1.3)),
                   week5_dabtram = log10(tab_mat_full[,"week5_DABTRAM"] + runif(n, min = 1, max = 1.3)),
                   chromvar1 = mean_mat[,j],
                   chromvar2 = max_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_mat)[j], "_mean")
  colnames(df)[5] <- paste0(colnames(max_mat)[j], "_max")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

###########

n <- nrow(tab_mat_full)
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v2-COCL2_ap1motif-day0size-day10size.pdf", 
    onefile = T, width = 8, height = 8)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat_full[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_cocl2 = log10(tab_mat_full[,"day10_COCL2"] + runif(n, min = 1, max = 1.3)),
                   week5_cocl2 = log10(tab_mat_full[,"week5_COCL2"] + runif(n, min = 1, max = 1.3)),
                   chromvar1 = mean_mat[,j],
                   chromvar2 = max_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_mat)[j], "_mean")
  colnames(df)[5] <- paste0(colnames(max_mat)[j], "_max")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

###########

n <- nrow(tab_mat_full)
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v2-CIS_ap1motif-day0size-day10size.pdf", 
    onefile = T, width = 8, height = 8)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat_full[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_cis = log10(tab_mat_full[,"day10_CIS"] + runif(n, min = 1, max = 1.3)),
                   week5_cis = log10(tab_mat_full[,"week5_CIS"] + runif(n, min = 1, max = 1.3)),
                   chromvar1 = mean_mat[,j],
                   chromvar2 = max_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_mat)[j], "_mean")
  colnames(df)[5] <- paste0(colnames(max_mat)[j], "_max")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

#############################

# check the PCA

set.seed(10)
n <- nrow(tab_mat_full)
x_vec <- log10(tab_mat_full[,"day10_CIS"] + runif(n, min = 1, max = 1.3))
all(rownames(max_mat) == rownames(tab_mat_full))
max_mat <- scale(max_mat)
pca_res <- stats::prcomp(max_mat)
motif_max_count <- pca_res$x[,1]

df <- data.frame(day10_cis = x_vec,
                 motif_max_pca = motif_max_count)
p <- GGally::ggpairs(df, 
                     lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                     progress = FALSE) + ggplot2::theme_bw()
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v2-CIS_ap1motif_pca_day10_day0motif.pdf", 
    onefile = T, width = 5, height = 5)
print(p)
dev.off()

###########################

# check the variance
tab_mat_short <- tab_mat_full[which(tab_mat_full[,"day0"] >= 3),]
lineage_idx_list2 <- lapply(rownames(tab_mat_short), function(lineage){
  intersect(which(all_data$assigned_lineage == lineage),
            which(all_data$dataset == "day0"))
})

sd_mat <- sapply(motif_idx, function(i){
  sapply(lineage_idx_list2, function(idx_vec){
    stats::sd(all_data[["chromvar"]]@data[i,idx_vec])
  })
})
rownames(sd_mat) <- rownames(tab_mat_short)
colnames(sd_mat) <- name_vec[motif_idx]

all(rownames(sd_mat) == rownames(tab_mat_short))

n <- nrow(tab_mat_short)
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v2-CIS_ap1motif_std-day0size-day10size.pdf", 
    onefile = T, width = 8, height = 8)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat_short[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_cis = log10(tab_mat_short[,"day10_CIS"] + runif(n, min = 1, max = 1.3)),
                   week5_cis = log10(tab_mat_short[,"week5_CIS"] + runif(n, min = 1, max = 1.3)),
                   chromvar = sd_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_mat)[j], "_std")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()
