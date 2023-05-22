rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_COCL2_spca.RData")

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec2 <- sort(unique(c(unlist(keygenes), keygenes_csc)))

gene_vec <- names(spca_res_list)
percent_rna <- sapply(gene_vec, function(gene){
  spca_res_list[[gene]]$U[1,1]^2
})
df <- data.frame(cv_score_vec = cv_score_vec,
                 gene = gene_vec,
                 labeling = gene_vec %in% gene_vec2,
                 percent_rna = percent_rna)
# put all the labeling == TRUE on bottom
df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = percent_rna, y = cv_score_vec))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                    ggplot2::aes(label = gene, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    max.overlaps = 50)
p1 <- p1 + Seurat::NoLegend()

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_keygenes-and-chrpeak_COCL2_rnapercentage_cvscore2.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

#################

gene_vec <- intersect(names(cv_score_vec)[which(cv_score_vec >= 0.44)],
                      names(percent_rna)[which(percent_rna <= 0.26)]) 
gene_vec <- sort(gene_vec)
length(gene_vec); length(gene_vec)/length(spca_res_list)

tmp_list <- lapply(gene_vec, function(gene){
  spca_res_list[[gene]]$dimred
})
spca_mat <- do.call(cbind, tmp_list)

##########

treatment <- "COCL2"
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
tmp <- metadata$assigned_lineage
names(tmp) <- rownames(metadata)
seurat_obj$assigned_lineage <- tmp[colnames(seurat_obj)]

set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:50) 

col_palette <- c("gray", "blue", "red")
names(col_palette) <- c("1loser_COCL2", "2mid_winner_COCL2", "3high_winner_COCL2")
names(col_palette) <- sort(unique(tier_vec))
p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                      cols = col_palette,
                      group.by = "tier_vec")
p1 <- p1 + ggplot2::ggtitle("COCL2 Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_COCL2_spca_day10_notrajectory2.png"),
                p1, device = "png", width = 7, height = 5, units = "in")


tab_mat <- tab_mat[order(tab_mat[,"week5_COCL2"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_COCL2"] >= 20)]

pdf("../../../../out/figures/Writeup6g/Writeup6g_COCL2_spca_day10_large-week5-lineages2.pdf", 
    onefile = T, width = 5, height = 5)
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  
  lineage_name <- lineage_vec[i]
  cell_idx <- intersect(which(metadata$assigned_lineage == lineage_name), 
                        which(metadata$dataset == "day10_COCL2"))
  cell_idx <- intersect(cell_idx,
                        which(metadata$assigned_posterior >= 0.5))
  cell_names <- rownames(metadata)[cell_idx]
  cell_names <- intersect(cell_names, rownames(rna_mat))
  
  p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                        cells.highlight = cell_names)
  p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: ", lineage_name, " (", length(cell_names), " day10 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()

##################################

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1.5)
alternative_clustering <- seurat_obj$seurat_clusters

set.seed(10)
seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1)

p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap", 
                      group.by = "RNA_snn_res.1.5", label = TRUE,
                      repel = TRUE, label.size = 2.5)
p1 <- p1 + ggplot2::ggtitle("COCL2 Day10 clusters")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_COCL2_spca_day10_clusters-finer.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap", 
                      group.by = "seurat_clusters", label = TRUE,
                      repel = TRUE, label.size = 2.5)
p1 <- p1 + ggplot2::ggtitle("COCL2 Day10 clusters")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_COCL2_spca_day10_clusters.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

tab_tier <- table(seurat_obj$seurat_clusters, seurat_obj$tier_vec)
tab_tier2 <- diag(1/Matrix::rowSums(tab_tier)) %*% tab_tier
rownames(tab_tier2) <- rownames(tab_tier)
round(100*tab_tier2)
tab_tier3 <- tab_tier %*% diag(1/Matrix::colSums(tab_tier))
colnames(tab_tier3) <- colnames(tab_tier)
round(100*tab_tier3)

##############################

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

treatment <- "COCL2"
relevant_celltypes <- c("day0", paste0(c("day10_", "week5_"), treatment))
all_data2 <- subset(all_data, dataset %in% relevant_celltypes)
Seurat::DefaultAssay(all_data2) <- "Saver"

all_data2[["geneActivity"]] <- NULL
all_data2[["fasttopic_CIS"]] <- NULL
all_data2[["fasttopic_DABTRAM"]] <- NULL
all_data2[["common_tcca"]] <- NULL
all_data2[["distinct1_tcca"]] <- NULL
all_data2[["distinct2_tcca"]] <- NULL

all_data2[["umap"]] <- NULL
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, 
                             reduction = "fasttopic_COCL2",
                             dims = 1:30)

# let's look into specifically Lin138675
lineage_name <- "Lin138675"
cell_names <- colnames(seurat_obj)[which(seurat_obj$assigned_lineage == lineage_name)]
all_data2[["umap"]]@cell.embeddings[cell_names,]
# the cell we want is day10_COCL2_TAGCCTTGTGTCCAGG-1
seurat_obj[["umap"]]@cell.embeddings[cell_names,]
seurat_obj[["umap"]]@cell.embeddings["day10_COCL2_TAGCCTTGTGTCCAGG-1",]
all_data2[["umap"]]@cell.embeddings["day10_COCL2_CGATTTGCAGTAGGAT-1",] 


# let's look into specifically Lin93927
lineage_name <- "Lin93927"
cell_names <- colnames(seurat_obj)[which(seurat_obj$assigned_lineage == lineage_name)]
all_data2[["umap"]]@cell.embeddings[cell_names,]
seurat_obj[["umap"]]@cell.embeddings[cell_names,]
all_data2[["umap"]]@cell.embeddings["day10_COCL2_CCAAATCAGAACCTGT-1",]
seurat_obj[["umap"]]@cell.embeddings["day10_COCL2_AATTGACGTAAGTCGC-1",]
seurat_obj[["umap"]]@cell.embeddings["day10_COCL2_TTATGCGCATAGCTTG-1",]

# plot a side-by-side of the various seurat_obj clusters
col_vec <- scales::hue_pal()(13)
names(col_vec) <- as.character(c(0:12))
pdf("../../../../out/figures/Writeup6g/Writeup6g_COCL2_spca_day10_clusters_side-by-side.pdf", 
    onefile = T, width = 10, height = 5)
for(cluster_num in 0:12){
  print(cluster_num)
  cell_names <- colnames(seurat_obj)[which(seurat_obj$seurat_clusters == cluster_num)]
  
  p1 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                        cells.highlight = cell_names,
                        cols.highlight = col_vec[as.character(cluster_num)])
  p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: All cells (Cluster ", cluster_num, ")\n",
                                     "Total #: ", length(cell_names), ", Loser%=", round(100*tab_tier2[as.character(cluster_num), "1loser_COCL2"]), 
                                     ", WinnerMid%=", round(100*tab_tier2[as.character(cluster_num), "2mid_winner_COCL2"]), 
                                     ", WinnerHigh%=", round(100*tab_tier2[as.character(cluster_num), "3high_winner_COCL2"])))
  p1 <- p1 + ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  
  p2 <- Seurat::DimPlot(seurat_obj, reduction = "umap", 
                        cells.highlight = cell_names,
                        cols.highlight = col_vec[as.character(cluster_num)])
  p2 <- p2 + ggplot2::ggtitle(paste0("COCL2: Day10 cells SPCA (Cluster ", cluster_num, ")\n",
                                     "%amongLoser=", round(100*tab_tier3[as.character(cluster_num), "1loser_COCL2"]), 
                                     ", %amongWinnerMid=", round(100*tab_tier3[as.character(cluster_num), "2mid_winner_COCL2"]), 
                                     ", %amongWinnerHigh=", round(100*tab_tier3[as.character(cluster_num), "3high_winner_COCL2"])))
  p2 <- p2 + ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  
  p <- cowplot::plot_grid(p1, p2, ncol = 2)
  print(p)
}
dev.off()