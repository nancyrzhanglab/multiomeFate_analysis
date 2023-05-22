rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)


load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_spca2.RData")
load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

# do the pseudotime scoring
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

# let's define winners and losers
# score each lineage by its mean day10 score
lineage_score <- sapply(1:nrow(tab_mat), function(i){
  if(i %% floor(nrow(tab_mat)/10) == 0) cat('*')
  lineage_name <- rownames(tab_mat)[i]
  cell_names <- colnames(seurat_obj)[which(seurat_obj$assigned_lineage == lineage_name)]
  if(length(cell_names) == 0) return(NA)
  mean(seurat_obj$pseudotime[cell_names], na.rm = T)
})
names(lineage_score) <- rownames(tab_mat)
table(is.na(lineage_score))
quantile(lineage_score, na.rm = T)
length(which(lineage_score >= 0.6))
length(which(lineage_score <= 0.3))

cell_winning_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score > 0.5)]),
                              which(all_data$dataset == "day0"))
cell_winning_idx <- intersect(cell_winning_idx,
                              which(all_data$assigned_posterior >= 0.5))
cell_losing_idx <- intersect(which(all_data$assigned_lineage %in% names(lineage_score)[which(lineage_score < 0.3)]),
                             which(all_data$dataset == "day0"))
cell_losing_idx <- intersect(cell_losing_idx,
                             which(all_data$assigned_posterior >= 0.5))
length(cell_winning_idx); length(cell_losing_idx)
keep_vec <- rep(NA, ncol(all_data))
keep_vec[cell_winning_idx] <- "winning_day0"
keep_vec[cell_losing_idx] <- "losing_day0"

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

motifs <- rownames(de_res)
data.use <- Signac::GetMotifData(object = all_data,
                                 assay = "ATAC",
                                 slot = "pwm")
data.use <- data.use[motifs]
names(data.use) <- Signac::GetMotifData(
  object = all_data,
  assay = "ATAC",
  slot = "motif.names"
)[motifs]
tmp <- de_res; rownames(tmp) <- names(data.use); tmp[1:50,]

names(data.use) <- sapply(1:length(data.use), function(i){
  paste0(names(data.use)[i], ", ", rownames(de_res)[i], 
         "\n-Log10pval=", round(-log10(de_res[i,"p_val"]), 2),
         ", value=", round(de_res[i,"avg_diff"], 2))
})
data.use <- data.use[1:28]

plot1 <- ggseqlogo::ggseqlogo(data = data.use, ncol = 4)
plot1 <- plot1 + ggplot2::theme_bw()
width <- 15; height <- 10

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_chromVar_day0_withday10-scores.png"),
                plot1, device = "png", width = width, height = height, units = "in")
