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

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
sum_vec <- tab_mat[,"day0"] + tab_mat[,"day10_COCL2"] + tab_mat[,"week5_COCL2"] 
tab_mat <- tab_mat[which(sum_vec >= 1), ]

cell_idx_list <- list(
  which(all_data$dataset == "day0"),
  which(all_data$dataset == "day10_COCL2"),
  which(all_data$dataset == "week5_COCL2")
)
names(cell_idx_list) <- c("day0", "day10_COCL2", "week5_COCL2")

# compute which cells belong to which lineage
lineage_names <- rownames(tab_mat)
lineage_idx_list_list <-
  lapply(cell_idx_list, function(idx_vec) {
    tmp <- lapply(lineage_names, function(lineage) {
      intersect(which(all_data$assigned_lineage == lineage),
                idx_vec)
    })
    names(tmp) <- lineage_names
    tmp
  })
names(lineage_idx_list_list) <- names(cell_idx_list)

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
      if(length(idx_vec) == 0) return(NA)
      mean(all_data[["chromvar"]]@data[i,idx_vec])
    })
  })
  
  rownames(mean_mat) <- names(lineage_idx_list)
  colnames(mean_mat) <- name_vec[motif_idx]
  
  mean_mat
})
names(mean_mat_list) <- names(lineage_idx_list_list)
all(rownames(mean_mat_list[[1]]) == rownames(mean_mat_list[[2]]))
all(colnames(mean_mat_list[[1]]) == colnames(mean_mat_list[[2]]))
all(rownames(mean_mat_list[[1]]) == rownames(tab_mat))

mean_all_mat <- sapply(motif_idx, function(i){
  sapply(1:length(lineage_idx_list_list[[1]]), function(j){
    idx_vec <- unlist(lapply(lineage_idx_list_list, function(lineage_idx_list){
      lineage_idx_list[[j]]
    }))
    stopifnot(length(unique(all_data$assigned_lineage[idx_vec])) == 1)
    mean(all_data[["chromvar"]]@data[i,idx_vec])
  })
})
rownames(mean_all_mat) <- lineage_names
colnames(mean_all_mat) <- name_vec[motif_idx]
all(rownames(mean_all_mat) == rownames(tab_mat))
all(colnames(mean_all_mat) == colnames(mean_mat_list[[1]]))

# make a ton of plots now
set.seed(10)
pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar-COCL2_ap1_each-timepoint-separate.pdf", 
    onefile = T, width = 10, height = 10)
n <- nrow(mean_mat_list[[1]])
p <- ncol(mean_mat_list[[1]])
for(j in 1:p){
  df <- data.frame(day0 = log10(tab_mat[,"day0"] + runif(n, min = 1, max = 1.3)),
                   day10_cocl2 = log10(tab_mat[,"day10_COCL2"] + runif(n, min = 1, max = 1.3)),
                   week5_cocl2 = log10(tab_mat[,"week5_COCL2"] + runif(n, min = 1, max = 1.3)),
                   chromvar1 = mean_mat_list[[1]][,j],
                   chromvar2 = mean_mat_list[[2]][,j],
                   chromvar3 = mean_mat_list[[3]][,j],
                   chromvar4 = mean_all_mat[,j])
  colnames(df)[4] <- paste0(colnames(mean_all_mat)[j], "_d0")
  colnames(df)[5] <- paste0(colnames(mean_all_mat)[j], "_d10")
  colnames(df)[6] <- paste0(colnames(mean_all_mat)[j], "_w5")
  colnames(df)[7] <- paste0(colnames(mean_all_mat)[j], "_all")
  p <- GGally::ggpairs(df, 
                       lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                       progress = FALSE) + ggplot2::theme_bw()
  print(p)
}
dev.off()
