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
  spca_res_list[[gene]]$dimred
})
spca_mat <- do.call(cbind, tmp_list)

##########

treatment <- "DABTRAM"
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
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
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)

p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap", 
                      group.by = "seurat_clusters", label = TRUE,
                      repel = TRUE, label.size = 2.5)
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 clusters")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_clusters.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

library(TrajectoryUtils)
# rm_clusters <- c("1", "3")
# keep_idx <- which(!seurat_obj$seurat_clusters %in% rm_clusters)
keep_idx <- 1:ncol(seurat_obj)

dist_method_vec = c("simple", "scaled.full", "scaled.diag", "slingshot", "mnn")

for(dist_method in dist_method_vec){
  mat <- seurat_obj[["pca"]]@cell.embeddings[keep_idx,1:30]
  mat <- scale(mat)
  print(dist_method)
  traj_res <- TrajectoryUtils:::.create_cluster_mst(
    x = mat,
    clusters = droplevels(factor(seurat_obj$seurat_clusters[keep_idx])),
    dist.method = dist_method
  )
  print(traj_res)
  print("=====")
}

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.5)
set.seed(10)
mat <- seurat_obj[["pca"]]@cell.embeddings[keep_idx,1:30]
mat <- scale(mat)
print(dist_method)
traj_res <- TrajectoryUtils:::.create_cluster_mst(
  x = mat,
  clusters = droplevels(factor(seurat_obj$seurat_clusters[keep_idx])),
  dist.method = "simple"
)
ordering_res <- TrajectoryUtils:::defineMSTPaths(
  g = traj_res,
  roots = "1"
)
ordering_res # this works, but is not great...

