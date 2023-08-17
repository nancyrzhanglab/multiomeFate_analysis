rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

# let's define winners and losers
# score each lineage by its mean day10 size
lineage_score <- tab_mat[,"day10_COCL2"]
names(lineage_score) <- rownames(tab_mat)

cell_winning_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score >= 10)]),
                              which(all_data$dataset == "day0"))
cell_winning_names <- colnames(all_data)[cell_winning_idx]
cell_losing_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score == 0)]),
                             which(all_data$dataset == "day0"))
cell_losing_names <- colnames(all_data)[cell_losing_idx]
keep_vec <- rep(NA, ncol(all_data))
names(keep_vec) <- colnames(all_data)
keep_vec[cell_winning_names] <- "winning_day0"
keep_vec[cell_losing_names] <- "losing_day0"

all_data$ident <- keep_vec
Seurat::Idents(all_data) <- "ident"

Seurat::DefaultAssay(all_data) <- "chromvar"
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = "winning_day0",
  ident.2 = "losing_day0",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  verbose = F
)

rowname_old <- rownames(de_res)
rowname_vec <- sapply(rownames(de_res), function(motif_name){
  idx <- which(rownames(all_data[["chromvar"]]@data) == motif_name)
  names(data.use)[idx]
})
rownames(de_res) <- rowname_vec

idx1 <- which(all_data$ident == "winning_day0")
idx2 <- which(all_data$ident == "losing_day0")
ttest_list <- lapply(1:nrow(all_data[["chromvar"]]@data), function(j){
  stats::t.test(
    x = all_data[["chromvar"]]@data[j, idx1],
    y = all_data[["chromvar"]]@data[j, idx2]
  )
})
names(ttest_list) <- names(data.use)

name_vec <- names(data.use) 
motif_idx <- sort(unique(c(grep("JUN", name_vec),
                           grep("FOS", name_vec),
                           grep("NFE2", name_vec),
                           grep("TEAD", name_vec))))
name_vec[motif_idx]

ttest_list[motif_idx[9]]

j <- motif_idx[9]
mean(all_data[["chromvar"]]@data[j, idx1])
mean(all_data[["chromvar"]]@data[j, idx2])

#############################

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- tab_mat[tab_mat[,"day0"] > 0,]
lineage_names <- rownames(tab_mat)

# compute which cells belong to which lineage
lineage_idx_list <- lapply(lineage_names, function(lineage){
  which(all_data$assigned_lineage == lineage)
})

name_vec <- names(data.use) 
mean_mat <- sapply(lineage_idx_list, function(idx_vec){
  Matrix::rowMeans(all_data[["chromvar"]]@data[,idx_vec,drop = F])
})
colnames(mean_mat) <- lineage_names
rownames(mean_mat) <- name_vec

lineage_winning <- rownames(tab_mat)[which(tab_mat[,"day10_COCL2"] >= 10)]
lineage_losing <- rownames(tab_mat)[which(tab_mat[,"day10_COCL2"] == 0)]

ttest_list <- lapply(1:nrow(mean_mat), function(j){
  stats::t.test(
    x = mean_mat[j, lineage_winning],
    y = mean_mat[j, lineage_losing]
  )
})
names(ttest_list) <- names(data.use)

ttest_list[motif_idx[9]]

# form a table
data_mat <- sapply(ttest_list, function(x){
  c(diff = x$estimate[1] - x$estimate[2], pvalue = x$p.value)
})
colnames(data_mat) <- names(data.use)
data_mat <- t(data_mat)
data_mat <- data_mat[order(data_mat[,"pvalue"], decreasing = F),]
data_mat[1:100,]

##############################

# let's try this. Let's make a scatter plot of the individual cell's chromvar against the average,
# and then we'll color points based on the status at day10

cells <- colnames(all_data)[intersect(
  which(all_data$dataset == "day0"),
  which(all_data$assigned_lineage %in% c(lineage_winning, lineage_losing))
)]

class_vec <- sapply(cells, function(cell_name){
  ifelse(all_data$assigned_lineage[cell_name] %in% lineage_winning, "win", "lose")
})
cell_chromvar <- as.numeric(all_data[["chromvar"]][which(name_vec == "TEAD3"), cells])
# compute the lineage chromvar from scratch just to double-check. This code is purposely slow
lineage_chromvar <- sapply(cells, function(cell_name){
  idx <- which(colnames(all_data[["chromvar"]]) == cell_name)
  lineage_name <- all_data$assigned_lineage[idx]
  cell_idx <- intersect(which(all_data$assigned_lineage == lineage_name),
                        which(all_data$dataset == "day0"))
  mean(all_data[["chromvar"]]@data[which(name_vec == "TEAD3"), cell_idx])
})
df <- data.frame(
  cell_chromvar = cell_chromvar,
  lineage_chromvar = lineage_chromvar,
  lineage_name = all_data$assigned_lineage[cells],
  class = class_vec
)
rownames(df) <- cells

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = cell_chromvar, y = lineage_chromvar, color = class_vec)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", linewidth = 1) + 
  ggplot2::labs(x = "Cell-wise chromVar", y = "Lineage-mean chromVar", title = "For TEAD3 (Stratified by COCL2 Day10 lineage size)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6l/Writeup6l_COCL2_TEAD3_day0-chromvar_by-day10-size.png"),
                p1, device = "png", width = 6, height = 5, units = "in")

#########

# I'm very confused. Isolate out the points that are in the "winning" lineages where the cell is above the average. 
# There should be one for each lineage, if the lineage is larger than 1
df2 <- df
df2 <- df2[df2$class == "win",]
df2 <- df2[df2$cell_chromvar > df2$lineage_chromvar + 1e-6,]
df2

tab_mat2 <- tab_mat[lineage_winning,]
tab_mat2 <- tab_mat2[order(tab_mat2[,"day0"], decreasing = T),]
lineage_name <- rownames(tab_mat2)[1]

cells <- colnames(all_data)[intersect(which(all_data$assigned_lineage == lineage_name),
                                      which(all_data$dataset == "day0"))]
df[cells,]
mean(df[cells,"cell_chromvar"])



