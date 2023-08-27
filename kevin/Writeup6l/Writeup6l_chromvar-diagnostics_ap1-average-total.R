rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(ggplot2)                     
library(GGally)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

stopifnot(length(data.use) == nrow(all_data[["chromvar"]]@data))
head(rownames(all_data[["chromvar"]]@data))
head(names(data.use))

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[tab_mat[,"day0"] > 0,]
lineage_names <- rownames(tab_mat)

name_vec <- names(data.use) 
motif_idx <- sort(unique(c(grep("JUN", name_vec),
                           grep("FOS", name_vec),
                           grep("NFE2", name_vec),
                           grep("TEAD", name_vec))))

# compute which cells belong to which lineage
lineage_idx_list <- lapply(lineage_names, function(lineage){
  which(all_data$assigned_lineage == lineage)
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
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar-avg-alltimepoints_DABTRAM_ap1motif-all-sizes.pdf", 
    onefile = T, width = 5, height = 5)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat[,"day0"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   day10_dabtram = log10(tab_mat[,"day10_DABTRAM"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   week5_dabtram = log10(tab_mat[,"week5_DABTRAM"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   chromvar = mean_mat[,j])
  colnames(df)[4] <- colnames(mean_mat)[j]
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

# make a ton of plots now
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar-avg-alltimepoints_CIS_ap1motif-all-sizes.pdf", 
    onefile = T, width = 5, height = 5)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat[,"day0"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   day10_cis = log10(tab_mat[,"day10_CIS"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   week5_cis = log10(tab_mat[,"week5_CIS"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   chromvar = mean_mat[,j])
  colnames(df)[4] <- colnames(mean_mat)[j]
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar-avg-alltimepoints_COCL2_ap1motif-all-sizes.pdf", 
    onefile = T, width = 5, height = 5)
for(j in 1:ncol(mean_mat)){
  df <- data.frame(day0 = log10(tab_mat[,"day0"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   day10_cocl2 = log10(tab_mat[,"day10_COCL2"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   week5_cocl2 = log10(tab_mat[,"week5_COCL2"] + runif(nrow(tab_mat), min = 1, max = 1.3)),
                   chromvar1 = mean_mat[,j])
  colnames(df)[4] <- colnames(mean_mat)[j]
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE)
  p <- p + ggplot2::theme_bw()
  print(p)
}
dev.off()
