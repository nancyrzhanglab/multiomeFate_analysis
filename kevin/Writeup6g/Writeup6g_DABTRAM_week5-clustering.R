rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activityPCA"]] <- NULL

treatment <- "DABTRAM"
all_data2 <- subset(all_data, dataset == paste0("week5_", treatment))
relevant_celltypes <- c("day0", paste0(c("day10_", "week5_"), treatment))
all_data <- subset(all_data, dataset %in% relevant_celltypes)

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_DABTRAM",
                            dims = 1:30)

############

set.seed(10)
all_data2 <- Seurat::FindNeighbors(all_data2, dims = 1:30, 
                                   reduction = "fasttopic_DABTRAM")
names(all_data2@graphs)
set.seed(10)
resolution <- 0.05
all_data2 <- Seurat::FindClusters(all_data2, 
                                  graph.name = "RNA_snn",
                                  resolution = resolution)
cluster_vec <- rep(NA, ncol(all_data))
names(cluster_vec) <- colnames(all_data)
cluster_vec[colnames(all_data2)] <- all_data2$seurat_clusters

all_data$custom_cluster <- cluster_vec

p1 <- Seurat::DimPlot(all_data, reduction = "umap", 
                      group.by = "custom_cluster", label = TRUE,
                      repel = TRUE, label.size = 2.5)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Week5 clusters,\nResolution: ", resolution))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_week5-clusters_2clusters.png"),
                p1, device = "png", width = 5, height = 5, units = "in")

######

# plot the % of lineages in each cluster against the size of the lineage
tab_mat2 <- table(all_data2$assigned_lineage, all_data2$seurat_clusters)

x_vec <- rowSums(tab_mat2)
y_vec <- tab_mat2[,"0"]/(tab_mat2[,"0"]+tab_mat2[,"1"])

names(x_vec) <- rownames(tab_mat2)
names(y_vec) <- rownames(tab_mat2)
rm_names <- names(x_vec)[sort(unique(c(which(x_vec == 0),
                                       which(is.nan(y_vec)))))]
if(length(rm_names) > 0){
  x_vec <- x_vec[which(!names(x_vec) %in% rm_names)]
  y_vec <- y_vec[which(!names(y_vec) %in% rm_names)]
}

df <- data.frame(size = log10(x_vec+1),
                 percentage = y_vec,
                 name = names(x_vec))
label_vec <- rep(F, length(x_vec))
idx1 <- which(x_vec >= quantile(x_vec, prob = 0.95))
label_vec[idx1] <- T
df$labeling = label_vec
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = size, y = percentage))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = name, color = labeling),
                                    size = 2,
                                    max.overlaps = 5)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Week5 clusters,\nPercentage in each cluster vs. size")) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_week5-clusters_percentage-vs-size.png"),
                p1, device = "png", width = 5, height = 5, units = "in")


