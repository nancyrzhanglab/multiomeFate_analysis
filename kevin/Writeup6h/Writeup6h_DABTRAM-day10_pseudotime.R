rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

all_data <- subset(all_data, dataset %in% c("day0", "day10_DABTRAM", "week5_DABTRAM"))
all_data <- subset(all_data, assigned_posterior >= 0.5)
all_data2 <- subset(all_data, dataset == "day10_DABTRAM")

###################

# do a pseudotime analysis just on the Day10s
Seurat::DefaultAssay(all_data2) <- "Saver"
all_data2 <- Seurat::RunPCA(all_data2, 
                            verbose = F,
                            assay = "Saver")
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, dims = 1:50) 

set.seed(10)
all_data2 <- Seurat::FindNeighbors(all_data2, dims = 1:30)
all_data2 <- Seurat::FindClusters(all_data2, resolution = 1.5)

set.seed(10)
all_data2 <- Seurat::FindClusters(all_data2, resolution = 1)

set.seed(10)
all_data2 <- Seurat::FindClusters(all_data2, resolution = 0.5)

set.seed(10)
all_data2 <- Seurat::FindClusters(all_data2, resolution = 0.25)

pdf(paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime-clusters.pdf"),
    onefile = T, width = 8, height = 6)
for(clustering in c("Saver_snn_res.1.5", "Saver_snn_res.1", "Saver_snn_res.0.5", "Saver_snn_res.0.25")){
  p <- Seurat::DimPlot(all_data2, reduction = "umap", 
                       group.by = clustering, label = TRUE,
                       repel = TRUE, label.size = 2.5)
  print(p)
}
dev.off() 

treatment <- "DABTRAM"
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]

tier1_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier1_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier1_names <- colnames(all_data)[tier1_idx]
tier2_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier2_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier2_names <- colnames(all_data)[tier2_idx]
tier3_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier3_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier3_names <- colnames(all_data)[tier3_idx]
tier_vec <- rep(NA, ncol(all_data))
names(tier_vec) <- colnames(all_data)
tier_vec[tier3_names] <- paste0("3high_winner_", treatment)
tier_vec[tier2_names] <- paste0("2mid_winner_", treatment)
tier_vec[tier1_names] <- paste0("1loser_", treatment)
all_data2$tier_vec <- tier_vec[colnames(all_data2)]

col_palette <- c("gray", "blue", "red")
names(col_palette) <- c("1loser_DABTRAM", "2mid_winner_DABTRAM", "3high_winner_DABTRAM")
names(col_palette) <- sort(unique(tier_vec))
p1 <- Seurat::DimPlot(all_data2, reduction = "umap",
                      cols = col_palette,
                      group.by = "tier_vec")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime-lineage-expansion.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

#############################################

source("../Writeup6g/slingshot_funcs.R")
order_vec <- c("0", "1", "2")
cluster_vec <- as.factor(as.character(all_data2$Saver_snn_res.0.25))
idx <- which(cluster_vec %in% order_vec)
initial_res <- .initial_curve_fit(
  cluster_vec = cluster_vec[idx],
  dimred = all_data2[["pca"]]@cell.embeddings[idx,1:30],
  lineage_order = order_vec
)
pseudotime_res <- .extract_pseudotime(
  dimred = all_data2[["pca"]]@cell.embeddings[,1:30],
  initial_fit = initial_res,
  stretch = 2
)
pseudotime_vec <- pseudotime_res/max(pseudotime_res)

res <- rep(NA, ncol(all_data))
names(res) <- colnames(all_data)
res[names(pseudotime_vec)] <- pseudotime_vec

# visualize it
all_data$pseudotime <- res
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_DABTRAM",
                            dims = 1:30)

p1 <- Seurat::DimPlot(all_data, reduction = "umap", 
                      group.by = "dataset")
p1 <- p1 + ggplot2::ggtitle("DABTRAM (Only cells w/ barcode)")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM-day10_dataset_all-cells.png"),
                p1, device = "png", width = 7, height = 5, units = "in")


p1 <- Seurat::FeaturePlot(all_data, reduction = "umap", 
                          features = "pseudotime")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day10 Pseudotime")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime_all-cells.png"),
                p1, device = "png", width = 7, height = 5, units = "in")


save(all_data, all_data2,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")


