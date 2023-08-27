rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(GGally)
library(ggplot2)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")
tab_mat_full <- table(all_data$assigned_lineage, all_data$dataset)
load("../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2_part4_lightweight-noATAC.RData")

tab_vec <- table(all_data$assigned_lineage)
tab_vec <- tab_vec[tab_vec > 0]
lineage_names <- names(tab_vec)
tab_mat_full <- tab_mat_full[lineage_names,]

name_vec <- names(data.use) 
motif_idx <- unique(c(grep("JUN", name_vec),
                      grep("FOS", name_vec),
                      grep("NFE2", name_vec),
                      grep("TEAD", name_vec)))

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

# make a ton of plots now
n <- nrow(tab_mat_full)
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v3-DABTRAM_ap1motif-day0size-day10size.pdf", 
    onefile = T, width = 8, height = 8)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat_full[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_dabtram = log10(tab_mat_full[,"day10_DABTRAM"] + runif(n, min = 1, max = 1.3)),
                   week5_dabtram = log10(tab_mat_full[,"week5_DABTRAM"] + runif(n, min = 1, max = 1.3)),
                   chromvar1 = mean_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_mat)[j], "_mean-day0")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

###########

n <- nrow(tab_mat_full)
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v3-COCL2_ap1motif-day0size-day10size.pdf", 
    onefile = T, width = 8, height = 8)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat_full[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_cocl2 = log10(tab_mat_full[,"day10_COCL2"] + runif(n, min = 1, max = 1.3)),
                   week5_cocl2 = log10(tab_mat_full[,"week5_COCL2"] + runif(n, min = 1, max = 1.3)),
                   chromvar1 = mean_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_mat)[j], "_mean-day0")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

###########

n <- nrow(tab_mat_full)
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_v3-CIS_ap1motif-day0size-day10size.pdf", 
    onefile = T, width = 8, height = 8)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat_full[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_cis = log10(tab_mat_full[,"day10_CIS"] + runif(n, min = 1, max = 1.3)),
                   week5_cis = log10(tab_mat_full[,"week5_CIS"] + runif(n, min = 1, max = 1.3)),
                   chromvar1 = mean_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_mat)[j], "_mean-day0")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

