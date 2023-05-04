rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
treatment <- "CIS"
gene_vec <- c("A1CF", "A3GALT2", "AAK1", "ABCA12", "ABCA13", "ABCA3", 
              "ABCA4", "ABCB1", "ABCB4", "ABCC1", "ABCC4", "ABCC9", 
              "ABCG1", "ABCG2", "ABHD17C", "ABI2")

surviving_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
dying_lineages <- rownames(tab_mat)[which(apply(tab_mat,1,max)<=1)]
winning_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% surviving_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
dying_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% dying_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
winning_cells <- colnames(all_data)[winning_idx]
dying_cells <- colnames(all_data)[dying_idx]

keep_vec <- rep(NA, ncol(all_data))
keep_vec[winning_idx] <- paste0("day0_winner_", treatment)
keep_vec[dying_idx] <- paste0("day0_loser_", treatment)
keep_vec <- as.factor(keep_vec)
table(keep_vec)
all_data2 <- all_data
all_data2$keep <- keep_vec
all_data2 <- subset(all_data2, keep %in% c(paste0("day0_winner_", treatment),
                                           paste0("day0_loser_", treatment)))

Seurat::DefaultAssay(all_data2) <- "ATAC"
Seurat::Idents(all_data2) <- "keep"

pdf(paste0("../../../../out/figures/Writeup6g/Writeup6g_coverage_", treatment, "_5000bp_tmp.pdf"),
    onefile = T, width = 9, height = 4.5)
for(gene in gene_vec){
  plot1 <- Signac::CoveragePlot(
    object = all_data2,
    region = gene,
    features = gene,
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  
  print(plot1)
}

dev.off() 
