rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_CIS.RData")

plot1 <- Seurat::DimPlot(all_data_subset, reduction = "scVelo.umap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("CIS: scVelo UMAP"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_scvelo_umap-cis.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

Seurat::DefaultAssay(all_data_subset) <- "RNA"
all_data_subset <- Seurat::FindNeighbors(all_data_subset, 
                                         reduction = "scVelo.umap", dims = 1:2)
all_data_subset <- Seurat::FindClusters(all_data_subset, resolution = 0.1)

# intersect with the dataset vector
clustering <- as.numeric(as.factor(all_data_subset$RNA_snn_res.0.1))
dataset <- as.numeric(as.factor(all_data_subset$dataset))
final_cluster <- clustering*max(dataset)+dataset
all_data_subset$final_cluster <- factor(final_cluster)

plot1 <- Seurat::DimPlot(all_data_subset, reduction = "scVelo.umap",
                         group.by = "RNA_snn_res.0.1", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("CIS: scVelo UMAP and clustering"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_scvelo_umap_clustering-cis.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


plot1 <- Seurat::DimPlot(all_data_subset, reduction = "scVelo.umap",
                         group.by = "final_cluster", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("CIS: scVelo UMAP and clustering"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_scvelo_umap_finalclustering-cis.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# ending clusters
idx <- grep("week5", all_data_subset$dataset)
end_clusters <- as.character(unique(all_data_subset$final_cluster[idx]))
set.seed(10)
slingshot_lin <- slingshot::getLineages(data = all_data_subset[["scVelo.umap"]]@cell.embeddings, 
                                        clusterLabels = as.character(all_data_subset$final_cluster),
                                        start.clus = c("10"))
slingshot_lin@metadata$lineages
slingshot_curve <- slingshot::getCurves(slingshot_lin)

pseudotime_mat <- slingshot_curve@assays@data$pseudotime
apply(pseudotime_mat, 2, function(x){quantile(x[!is.na(x)])})
pseudotime_mat <- apply(pseudotime_mat, 2, function(x){
  idx <- which(!is.na(x))
  x[idx] <- x[idx]/max(x[idx])
  x
})
pseudotime_vec <- apply(pseudotime_mat, 1, function(x){
  x <- x[!is.na(x)]
  mean(x)
})
all_data_subset$pseudotime <- pseudotime_vec

plot1 <- Seurat::FeaturePlot(all_data_subset, reduction = "scVelo.umap",
                         features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("CIS: scVelo UMAP, Slingshot pseudotime"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_scvelo_umap_slingshot-pseudotime-cis.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

curves_list <- slingshot_curve@metadata$curves
png("../../../../out/figures/Writeup4d/Writeup4d_scvelo_umap_slingshot-curves-cis.png",
    height = 3000, width = 3000, units = "px", res = 300)
plot(x = all_data_subset[["scVelo.umap"]]@cell.embeddings[,1],
     y = all_data_subset[["scVelo.umap"]]@cell.embeddings[,2],
     pch = 16, col = as.numeric(as.factor(all_data_subset$dataset)))
for(i in 1:length(curves_list)){
  ord <- curves_list[[i]]$ord
  mat <- curves_list[[i]]$s[ord,]
  lines(mat[,1], mat[,2], lwd = 10, col = "white")
  lines(mat[,1], mat[,2], lwd = 7)
}
graphics.off()
