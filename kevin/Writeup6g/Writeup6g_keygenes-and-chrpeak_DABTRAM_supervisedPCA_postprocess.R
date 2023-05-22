rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_spca2.RData")

# make the embedding, and see if lineages coalesce in the embedding
gene_vec <- names(spca_res_list)
percent_rna <- sapply(gene_vec, function(gene){
  spca_res_list[[gene]]$U[1,1]^2
})

round(quantile(100*cv_score_vec, probs = seq(0,1,length.out=11)))
round(quantile(100*percent_rna, probs = seq(0,1,length.out=11)))

gene_vec <- intersect(names(cv_score_vec)[which(cv_score_vec >= 0.5)],
                      names(percent_rna)[which(percent_rna <= 0.3)]) 
gene_vec <- sort(gene_vec)
length(gene_vec); length(gene_vec)/length(spca_res_list)

tmp_list <- lapply(gene_vec, function(gene){
  Re(spca_res_list[[gene]]$dimred)
})
spca_mat <- do.call(cbind, tmp_list)

##########

treatment <- "DABTRAM"
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)

tier1_idx <- intersect(
  intersect(which(metadata$assigned_lineage %in% tier1_lineages),
            which(metadata$assigned_posterior >= 0.5)),
  which(metadata$dataset == paste0("day10_", treatment))
)
tier1_names <- intersect(rownames(metadata)[tier1_idx], rownames(rna_mat))
tier2_idx <- intersect(
  intersect(which(metadata$assigned_lineage %in% tier2_lineages),
            which(metadata$assigned_posterior >= 0.5)),
  which(metadata$dataset == paste0("day10_", treatment))
)
tier2_names <- intersect(rownames(metadata)[tier2_idx], rownames(rna_mat))
tier3_idx <- intersect(
  intersect(which(metadata$assigned_lineage %in% tier3_lineages),
            which(metadata$assigned_posterior >= 0.5)),
  which(metadata$dataset == paste0("day10_", treatment))
)
tier3_names <- intersect(rownames(metadata)[tier3_idx], rownames(rna_mat))
tier_vec <- rep(NA, nrow(rna_mat))
names(tier_vec) <- rownames(rna_mat)
tier_vec[tier3_names] <- paste0("3high_winner_", treatment)
tier_vec[tier2_names] <- paste0("2mid_winner_", treatment)
tier_vec[tier1_names] <- paste0("1loser_", treatment)

table(tier_vec); table(is.na(tier_vec))

##########

set.seed(10)
seurat_obj <- Seurat::CreateSeuratObject(counts = t(spca_mat))
seurat_obj[["RNA"]]@data <- seurat_obj[["RNA"]]@counts
seurat_obj[["RNA"]]@scale.data <- as.matrix(seurat_obj[["RNA"]]@counts)
seurat_obj[["RNA"]]@var.features <- rownames(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = F)
seurat_obj$tier_vec <- tier_vec
tmp <- metadata$assigned_lineage
names(tmp) <- rownames(metadata)
seurat_obj$assigned_lineage <- tmp[colnames(seurat_obj)]

set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:50) 

col_palette <- c("gray", "blue", "red")
names(col_palette) <- c("1loser_DABTRAM", "2mid_winner_DABTRAM", "3high_winner_DABTRAM")
names(col_palette) <- sort(unique(tier_vec))
p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                      cols = col_palette,
                      group.by = "tier_vec")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_notrajectory2.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

tab_mat <- tab_mat[order(tab_mat[,"week5_DABTRAM"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_DABTRAM"] >= 20)]

pdf("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_large-week5-lineages2.pdf", 
    onefile = T, width = 5, height = 5)
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  
  lineage_name <- lineage_vec[i]
  cell_idx <- intersect(which(metadata$assigned_lineage == lineage_name), 
                        which(metadata$dataset == "day10_DABTRAM"))
  cell_idx <- intersect(cell_idx,
                        which(metadata$assigned_posterior >= 0.5))
  cell_names <- rownames(metadata)[cell_idx]
  cell_names <- intersect(cell_names, rownames(rna_mat))
  
  p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                        cells.highlight = cell_names)
  p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: ", lineage_name, " (", length(cell_names), " day10 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()

###############

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1.5)
alternative_clustering <- seurat_obj$seurat_clusters

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1)

# fix specifically cluster 4
tmp <- table(seurat_obj$RNA_snn_res.1, seurat_obj$RNA_snn_res.1.5)
idx <- intersect(which(seurat_obj$RNA_snn_res.1 == "4"),
                 which(seurat_obj$RNA_snn_res.1.5 %in% c("0", "5", "8", "13")))
seurat_obj$seurat_clusters[idx] <- NA

p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap", 
                      group.by = "RNA_snn_res.1.5", label = TRUE,
                      repel = TRUE, label.size = 2.5)
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 clusters")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_clusters-finer.png"),
                p1, device = "png", width = 7, height = 5, units = "in")


p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap", 
                      group.by = "seurat_clusters", label = TRUE,
                      repel = TRUE, label.size = 2.5)
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 clusters")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_clusters.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

# let's try something else: compute the percentage of cells by their tiers
tab_tier <- table(seurat_obj$seurat_clusters, seurat_obj$tier_vec)
tab_tier <- diag(1/Matrix::rowSums(tab_tier)) %*% tab_tier
round(100*tab_tier)

# let's make the order manually for now
source("slingshot_funcs.R")
order_list <- list(c("2", "1", "4", "7", "6", "10"),
                   c("2", "9", "0", "3", "5", "8"))
cluster_vec <- as.factor(as.character(seurat_obj$seurat_clusters))
pseudotime_mat <- sapply(1:2, function(k){
  order_vec <- order_list[[k]]
  
  idx <- which(cluster_vec %in% order_vec)
  initial_res <- .initial_curve_fit(
    cluster_vec = cluster_vec[idx],
    dimred = seurat_obj[["pca"]]@cell.embeddings[idx,1:30],
    lineage_order = order_vec
  )
  
  pseudotime_res <- .extract_pseudotime(
    dimred = seurat_obj[["pca"]]@cell.embeddings[idx,1:30],
    initial_fit = initial_res,
    stretch = 2
  )
  pseudotime_res <- pseudotime_res/max(pseudotime_res)
  
  res <- rep(NA, ncol(seurat_obj))
  res[idx] <- pseudotime_res
  # res[idx] <- rank(pseudotime_res)/length(idx)
  res
})
pseudotime_vec <- apply(pseudotime_mat, 1, function(x){
  mean(x, na.rm=T)
})
# fill in the NAs
na_idx <- which(is.na(pseudotime_vec))
if(length(na_idx) > 0){
  graph <- seurat_obj@graphs[["RNA_snn"]]
  pseudotime_vec2 <- pseudotime_vec
  for(i in na_idx){
    neighbor_idx <- multiomeFate:::.nonzero_col(
      mat = graph, 
      col_idx = i, 
      bool_value = F
    )
    stopifnot(any(!is.na(pseudotime_vec[neighbor_idx])))
    pseudotime_vec2[i] <- mean(pseudotime_vec[neighbor_idx], na.rm = T)
  }
  
  pseudotime_vec <- pseudotime_vec2
}
stopifnot(all(!is.na(pseudotime_vec)))
seurat_obj$pseudotime <- pseudotime_vec

p1 <- Seurat::FeaturePlot(seurat_obj, reduction = "umap", 
                          features = "pseudotime")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 Pseudotime")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_pseudotime.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

####

mat <- seurat_obj[["umap"]]@cell.embeddings
pseudotime_vec <- seurat_obj$pseudotime
color_palette_vec <- grDevices::colorRampPalette(c("lightgray", "blue"))(100)
spacing_vec <- seq(0,1,length.out=100)
color_vec <- sapply(pseudotime_vec, function(val){
  color_palette_vec[which.min(abs(val - spacing_vec))]
})

png("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_pseudotime_split.png",
    height = 1500, width = 3000, units = "px", res = 300)
par(mfrow = c(1,2))
for(k in 1:2){
  cell_idx <- which(seurat_obj$seurat_clusters %in% order_list[[k]])
  plot(mat[,1], mat[,2], xlab = "UMAP 1", ylab = "UMAP 2", 
       main = paste0("Lineage ", k),
       pch = 16, col = "gray90", cex = 0.5)
  points(mat[cell_idx,1], mat[cell_idx,2], pch = 16, col = color_vec[cell_idx], cex = 1)
}
graphics.off()

########################################

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

treatment <- "DABTRAM"
relevant_celltypes <- c("day0", paste0(c("day10_", "week5_"), treatment))
all_data2 <- subset(all_data, dataset %in% relevant_celltypes)
Seurat::DefaultAssay(all_data2) <- "Saver"

all_data2[["geneActivity"]] <- NULL
all_data2[["fasttopic_CIS"]] <- NULL
all_data2[["fasttopic_COCL2"]] <- NULL
all_data2[["common_tcca"]] <- NULL
all_data2[["distinct1_tcca"]] <- NULL
all_data2[["distinct2_tcca"]] <- NULL
all_data2[["activityPCA"]] <- NULL

all_data2[["umap"]] <- NULL
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, 
                             reduction = "fasttopic_DABTRAM",
                             dims = 1:30)

###

# plot week5 cells by their lineage
tab_mat <- table(all_data2$assigned_lineage, all_data2$dataset)
length(which(tab_mat[,"week5_DABTRAM"] >= 20))
cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")

tab_mat <- tab_mat[order(tab_mat[,"week5_DABTRAM"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_DABTRAM"] >= 20)]

pseudotime_vec <- seurat_obj$pseudotime
color_pseudotime_palette_vec <- grDevices::colorRampPalette(c("lightgray", "blue"))(100)
pseudotime_spacing_vec <- seq(0,1,length.out=100)
color_pseudotime_vec <- sapply(pseudotime_vec, function(val){
  color_pseudotime_palette_vec[which.min(abs(val - pseudotime_spacing_vec))]
})
names(color_pseudotime_vec) <- colnames(seurat_obj)

tab_mat <- tab_mat[order(tab_mat[,"week5_DABTRAM"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_DABTRAM"] >= 20)]
pdf("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_pseudotime_large-lineages_week5.pdf", 
    onefile = T, width = 15, height = 5)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  lineage_name <- lineage_vec[i]
  cell_idx <- intersect(which(all_data2$dataset == paste0("week5_", treatment)),
                        which(all_data2$assigned_lineage == lineage_name))
  cell_idx <- intersect(cell_idx, 
                        which(all_data2$assigned_posterior >= 0.5))
  cell_names <- colnames(all_data2)[cell_idx]
  
  cell_idx2 <- intersect(which(all_data2$dataset == paste0("day10_", treatment)),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx2 <- intersect(cell_idx2, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names2 <- colnames(all_data2)[cell_idx2]
  
  cell_idx3 <- intersect(which(all_data2$dataset == "day0"),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx3 <- intersect(cell_idx3, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names3 <- colnames(all_data2)[cell_idx3]
  
  mat <- all_data2[["umap"]]@cell.embeddings
  plot(mat[,1], mat[,2], xlab = "UMAP 1", ylab = "UMAP 2", 
       main = paste0("DABTRAM: ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                     length(cell_idx2), " day10 cells and ", length(cell_idx3), " day0 cells)"),
       pch = 16, col = "gray90", cex = 0.5)
  points(mat[cell_names,1], mat[cell_names,2], pch = 16, col = cols.highlight_vec["week5"], cex = 1.5)
  points(mat[cell_names2,1], mat[cell_names2,2], pch = 16, col = cols.highlight_vec["day10"], cex = 1.5)
  points(mat[cell_names3,1], mat[cell_names3,2], pch = 16, col = cols.highlight_vec["day0"], cex = 1.5)
  
  ##
  
  mat <- seurat_obj[["umap"]]@cell.embeddings
  for(k in 1:2){
    if(k == 1) {
      cell_idx <- intersect(which(seurat_obj$seurat_clusters %in% order_list[[k]]),
                            which(seurat_obj$assigned_lineage == lineage_name))
    } else {
      cell_idx <- intersect(c(which(seurat_obj$seurat_clusters %in% order_list[[k]]),
                              which(is.na(seurat_obj$seurat_clusters))),
                            which(seurat_obj$assigned_lineage == lineage_name))
    }
    
    plot(mat[,1], mat[,2], xlab = "UMAP 1", ylab = "UMAP 2", 
         main = paste0("Day10 trajectory ", k),
         pch = 16, col = "gray90", cex = 0.5)
    points(mat[cell_idx,1], mat[cell_idx,2], pch = 16, col = color_pseudotime_vec[cell_idx], cex = 1.5)
  }
  
}
dev.off()

################

png("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_pseudotime_all-timepoints.png",
    height = 2500, width = 2500, units = "px", res = 300)

mat <- all_data2[["umap"]]@cell.embeddings
plot(mat[,1], mat[,2], xlab = "UMAP 1", ylab = "UMAP 2", 
     main = paste0("DABTRAM Day10 pseudotime\namong all timepoints"),
     pch = 16, col = "gray90", cex = 0.5)
cell_names <- colnames(seurat_obj)
points(mat[cell_names,1], mat[cell_names,2], pch = 16, col = color_pseudotime_vec[cell_names])
graphics.off()

df <- data.frame(tier = seurat_obj$tier_vec,
                 pseudotime = seurat_obj$pseudotime)
col_vec <- c("lightgray", "#7190AF", "dodgerblue4")
names(col_vec) <- c("1loser_DABTRAM", "2mid_winner_DABTRAM", "3high_winner_DABTRAM")

p1 <- ggplot2::ggplot(df, ggplot2::aes(x=tier, y=pseudotime)) +
  ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=tier)) +
  ggplot2::scale_fill_manual(values=col_vec) +
  ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::geom_boxplot(width=0.05)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_pseudotime_violinplot.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

##############

tab_tier <- table(seurat_obj$seurat_clusters, seurat_obj$tier_vec)
tab_tier2 <- diag(1/Matrix::rowSums(tab_tier)) %*% tab_tier
rownames(tab_tier2) <- rownames(tab_tier)
round(100*tab_tier2)
tab_tier3 <- tab_tier %*% diag(1/Matrix::colSums(tab_tier))
colnames(tab_tier3) <- colnames(tab_tier)
round(100*tab_tier3)

# plot a side-by-side of the various seurat_obj clusters
col_vec <- scales::hue_pal()(11)
names(col_vec) <- as.character(c(0:10))
pdf("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_clusters_side-by-side.pdf", 
    onefile = T, width = 10, height = 5)
for(cluster_num in 0:10){
  print(cluster_num)
  cell_names <- colnames(seurat_obj)[which(seurat_obj$seurat_clusters == cluster_num)]
  
  p1 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                        cells.highlight = cell_names,
                        cols.highlight = col_vec[as.character(cluster_num)])
  p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: All cells (Cluster ", cluster_num, ")\n",
                                     "Total #: ", length(cell_names), ", Loser%=", round(100*tab_tier2[as.character(cluster_num), "1loser_DABTRAM"]), 
                                     ", WinnerMid%=", round(100*tab_tier2[as.character(cluster_num), "2mid_winner_DABTRAM"]), 
                                     ", WinnerHigh%=", round(100*tab_tier2[as.character(cluster_num), "3high_winner_DABTRAM"])))
  p1 <- p1 + ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  
  p2 <- Seurat::DimPlot(seurat_obj, reduction = "umap", 
                        cells.highlight = cell_names,
                        cols.highlight = col_vec[as.character(cluster_num)])
  p2 <- p2 + ggplot2::ggtitle(paste0("DABTRAM: Day10 cells SPCA (Cluster ", cluster_num, ")\n",
                                     "%amongLoser=", round(100*tab_tier3[as.character(cluster_num), "1loser_DABTRAM"]), 
                                     ", %amongWinnerMid=", round(100*tab_tier3[as.character(cluster_num), "2mid_winner_DABTRAM"]), 
                                     ", %amongWinnerHigh=", round(100*tab_tier3[as.character(cluster_num), "3high_winner_DABTRAM"])))
  p2 <- p2 + ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  
  p <- cowplot::plot_grid(p1, p2, ncol = 2)
  print(p)
}
dev.off()

########################################

## NOTE: Junk code since I couldn't figure out using existing tools how to figure out the start and ends
# library(TrajectoryUtils)
# # rm_clusters <- c("1", "3")
# # keep_idx <- which(!seurat_obj$seurat_clusters %in% rm_clusters)
# keep_idx <- 1:ncol(seurat_obj)
# 
# dist_method_vec = c("simple", "scaled.full", "scaled.diag", "slingshot", "mnn")
# 
# for(dist_method in dist_method_vec){
#   mat <- seurat_obj[["pca"]]@cell.embeddings[keep_idx,1:30]
#   # mat <- scale(mat)
#   print(dist_method)
#   traj_res <- TrajectoryUtils:::.create_cluster_mst(
#     x = mat,
#     clusters = droplevels(factor(seurat_obj$seurat_clusters[keep_idx])),
#     dist.method = dist_method
#   )
#   print(traj_res)
#   print("=====")
# }
# 
# set.seed(10)
# seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
# seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1)
# set.seed(10)
# mat <- seurat_obj[["pca"]]@cell.embeddings[keep_idx,1:30]
# # mat <- scale(mat)
# traj_res <- TrajectoryUtils:::.create_cluster_mst(
#   x = mat,
#   clusters = droplevels(factor(seurat_obj$seurat_clusters[keep_idx])),
#   dist.method = "simple"
# )
# ordering_res <- TrajectoryUtils:::defineMSTPaths(
#   g = traj_res,
#   roots = "2"
# )
# ordering_res 
