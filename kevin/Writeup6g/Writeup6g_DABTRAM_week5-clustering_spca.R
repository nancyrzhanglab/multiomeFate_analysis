rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_spca2.RData")

####################################################################
# let's make the spca embedding first 
# make the embedding, and see if lineages coalesce in the embedding
gene_vec <- names(spca_res_list)
percent_rna <- sapply(gene_vec, function(gene){
  spca_res_list[[gene]]$U[1,1]^2
})

gene_vec <- intersect(names(cv_score_vec)[which(cv_score_vec >= 0.5)],
                      names(percent_rna)[which(percent_rna <= 0.3)]) 
gene_vec <- sort(gene_vec)

tmp_list <- lapply(gene_vec, function(gene){
  Re(spca_res_list[[gene]]$dimred)
})
spca_mat <- do.call(cbind, tmp_list)

treatment <- "DABTRAM"
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]

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

####################################################################

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL

treatment <- "DABTRAM"
all_data2 <- subset(all_data, dataset == paste0("week5_", treatment))
relevant_celltypes <- c("day0", paste0(c("day10_", "week5_"), treatment))
all_data <- subset(all_data, dataset %in% relevant_celltypes)

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

###########

tab_mat2 <- table(all_data2$assigned_lineage, all_data2$seurat_clusters)
tab_mat3 <- table(seurat_obj$assigned_lineage, seurat_obj$seurat_clusters)

tab_mat2 <- tab_mat2[,c("0", "1")]
tab_mat2 <- tab_mat2[rowSums(tab_mat2) > 1, ]
tab_mat3 <- tab_mat3[rowSums(tab_mat3) > 1, ]
large_lineages <- rownames(tab_mat2)[which(rowSums(tab_mat2) >= 20)]

lineage_day10_percentage <- sapply(1:nrow(tab_mat3), function(i){
  x <- tab_mat3[i,]
  percentage_vec <- sapply(order_list, function(order_vec){
    sum(x[order_vec])
  })
  percentage_vec/sum(percentage_vec)
})
lineage_day10_percentage <- t(lineage_day10_percentage)
rownames(lineage_day10_percentage) <- rownames(tab_mat3)

shared_lineages <- sort(intersect(rownames(lineage_day10_percentage), rownames(tab_mat2)))
shared_lineages <- intersect(shared_lineages, large_lineages)
length(shared_lineages)

lineage_day10_percentage <- lineage_day10_percentage[shared_lineages,]
lineage_week5_percentage <- tab_mat2[shared_lineages,]
for(i in 1:nrow(lineage_week5_percentage)){
  lineage_week5_percentage[i,] <- lineage_week5_percentage[i,]/sum(lineage_week5_percentage[i,])
}

# figure out which is the pairing of clusters
cor_1 <- stats::cor(lineage_day10_percentage[,1], lineage_week5_percentage[,1])
cor_2 <- stats::cor(lineage_day10_percentage[,1], lineage_week5_percentage[,2])
cor_1; cor_2

if(cor_2 > cor_1){
  lineage_week5_percentage <- lineage_week5_percentage[,c(2,1)]
}

df <- data.frame(day10_percentage = lineage_day10_percentage[,1],
                 week5_percentage = lineage_week5_percentage[,1],
                 name = rownames(lineage_day10_percentage),
                 count = sqrt(rowSums(tab_mat3)[rownames(lineage_day10_percentage)]))

p1 <- ggplot2::ggplot(df, ggplot2::aes(x = day10_percentage, y = week5_percentage))
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(size = count))
p1 <- p1 + ggrepel::geom_text_repel(data = df,
                                    ggplot2::aes(label = name),
                                    size = 2,
                                    max.overlaps = 5)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM Day10-Week5 percentage\nCorrelation: ", round(max(cor_1,cor_2),2))) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_day10-week5_percentage.png"),
                p1, device = "png", width = 5, height = 5, units = "in")




