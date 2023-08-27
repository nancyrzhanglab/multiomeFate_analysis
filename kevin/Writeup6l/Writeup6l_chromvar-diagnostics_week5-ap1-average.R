rm(list = ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(GGally)

load(
  "../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData"
)

# see if the AP1 motifs at week5 are correlated with their lineage size
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

cell_week5_idx_list <- list(
  which(all_data$dataset == "week5_CIS"),
  which(all_data$dataset == "week5_COCL2"),
  which(all_data$dataset == "week5_DABTRAM")
)
names(cell_week5_idx_list) <- c("CIS", "COCL2", "DABTRAM")

# compute which cells belong to which lineage
lineage_names <- rownames(tab_mat)
lineage_idx_list_list <-
  lapply(cell_week5_idx_list, function(idx_vec) {
    tmp <- lapply(lineage_names, function(lineage) {
      intersect(which(all_data$assigned_lineage == lineage),
                idx_vec)
    })
    names(tmp) <- lineage_names
    tmp[sapply(tmp, length) > 0]
  })
names(lineage_idx_list_list) <- names(cell_week5_idx_list)
# just a test
length(unique(all_data$dataset[lineage_idx_list_list[[1]][[43]]])) == 1
length(unique(all_data$assigned_lineage[lineage_idx_list_list[[1]][[43]]])) == 1

name_vec <- names(data.use)
motif_idx <- unique(c(
  grep("JUN", name_vec),
  grep("FOS", name_vec),
  grep("NFE2", name_vec),
  grep("TEAD", name_vec)
))

# for each lineage, compute the mean scores
mean_mat_list <- lapply(lineage_idx_list_list, function(lineage_idx_list){
  mean_mat <- sapply(motif_idx, function(i){
    sapply(lineage_idx_list, function(idx_vec){
      mean(all_data[["chromvar"]]@data[i,idx_vec])
    })
  })
  
  rownames(mean_mat) <- names(lineage_idx_list)
  colnames(mean_mat) <- name_vec[motif_idx]
  
  mean_mat
})
names(mean_mat_list) <- names(lineage_idx_list_list)

# make a ton of plots now
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar-CIS_ap1motif-byweek5_week5size.pdf", 
    onefile = T, width = 5, height = 5)
mean_mat <- mean_mat_list[["CIS"]]
lineage_name_tmp <- rownames(mean_mat)
p <- ncol(mean_mat)
for(j in 1:p){
  df <- data.frame(day0 = log10(tab_mat[lineage_name_tmp,"day0"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   day10_cis = log10(tab_mat[lineage_name_tmp,"day10_CIS"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   week5_cis = log10(tab_mat[lineage_name_tmp,"week5_CIS"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   chromvar = mean_mat[,j])
  colnames(df)[4] <- colnames(mean_mat)[j]
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar-COCL2_ap1motif-byweek5_week5size.pdf", 
    onefile = T, width = 5, height = 5)
mean_mat <- mean_mat_list[["COCL2"]]
lineage_name_tmp <- rownames(mean_mat)
p <- ncol(mean_mat)
for(j in 1:p){
  df <- data.frame(day0 = log10(tab_mat[lineage_name_tmp,"day0"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   day10_cocl2 = log10(tab_mat[lineage_name_tmp,"day10_COCL2"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   week5_cocl2 = log10(tab_mat[lineage_name_tmp,"week5_COCL2"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   chromvar = mean_mat[,j])
  colnames(df)[4] <- colnames(mean_mat)[j]
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()

set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar-DABTRAM_ap1motif-byweek5_week5size.pdf", 
    onefile = T, width = 5, height = 5)
mean_mat <- mean_mat_list[["DABTRAM"]]
lineage_name_tmp <- rownames(mean_mat)
p <- ncol(mean_mat)
for(j in 1:p){
  df <- data.frame(day0 = log10(tab_mat[lineage_name_tmp,"day0"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   day10_dabtram = log10(tab_mat[lineage_name_tmp,"day10_DABTRAM"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   week5_dabtram = log10(tab_mat[lineage_name_tmp,"week5_DABTRAM"] + runif(length(lineage_name_tmp), min = 1, max = 1.3)),
                   chromvar = mean_mat[,j])
  colnames(df)[4] <- colnames(mean_mat)[j]
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()



