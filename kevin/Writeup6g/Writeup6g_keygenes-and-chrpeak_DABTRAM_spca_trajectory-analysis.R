rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)
library(ordinal)

load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_spca.RData")

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

set.seed(10)
seurat_obj <- Seurat::CreateSeuratObject(counts = t(spca_mat))
seurat_obj[["RNA"]]@data <- seurat_obj[["RNA"]]@counts
seurat_obj[["RNA"]]@scale.data <- as.matrix(seurat_obj[["RNA"]]@counts)
seurat_obj[["RNA"]]@var.features <- rownames(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = F)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:50)
seurat_obj$tier_vec <- tier_vec

col_palette <- c("gray", "blue", "red")
names(col_palette) <- c("1loser_DABTRAM", "2mid_winner_DABTRAM", "3high_winner_DABTRAM")
names(col_palette) <- sort(unique(tier_vec))
p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                      cols = col_palette,
                      group.by = "tier_vec")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_notrajectory.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

tab_mat <- tab_mat[order(tab_mat[,"week5_DABTRAM"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_DABTRAM"] >= 20)]

pdf("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10_large-week5-lineages.pdf", 
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
  p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: ", lineage_name, " (", length(cell_names), " day10 cells")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()

