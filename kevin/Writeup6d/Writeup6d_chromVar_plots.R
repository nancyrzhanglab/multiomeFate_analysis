rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")
source("../Writeup6b/gene_list.R")

Seurat::DefaultAssay(all_data) <- "chromvar"

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
gene_vec <- sort(unique(unlist(keygenes)))
gene_vec <- gene_vec[which(gene_vec %in% rownames(all_data[["RNA"]]))]

###########################################
# let's just focus on DABTRAM
# first, find all those differential motifs at day0

treatment <- "DABTRAM"
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
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
length(winning_idx); length(dying_idx)
ident_vec <- rep(NA, ncol(all_data))
ident_vec[winning_idx] <- paste0("day0_win_", treatment)
ident_vec[dying_idx] <- paste0("day0_lose_", treatment)
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"
table(Seurat::Idents(all_data))
all_data <- subset(all_data, ident %in% c(paste0("day0_win_", treatment), paste0("day0_lose_", treatment)))

# DE analysis
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = paste0("day0_win_", treatment),
  ident.2 = paste0("day0_lose_", treatment),
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  verbose = F
)

idx <- which(de_res[,"p_val_adj"] <= 1e-2)
de_res[idx,]

# get the motif names
motifs <- rownames(de_res)[idx]
motif_names <- unlist(Signac::GetMotifData(
  object = all_data,
  assay = "ATAC",
  slot = "motif.names"
)[motifs])
motif_names

# JANKY: Do surgery and change the rownames
stopifnot(!any(duplicated(motif_names)))
rowname_vec <- rownames(all_data[["chromvar"]]@data)
for(i in 1:length(motifs)){
  rowname_vec[which(rowname_vec == motifs[i])] <- motif_names[i]
}
rownames(all_data[["chromvar"]]@data) <- rowname_vec

# make some violin plots
plot1 <- Seurat::VlnPlot(all_data, 
                         features = motif_names,
                         slot = "data",
                         ncol = 6,
                         pt.size = 0.5)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6d/Writeup6d_DABTRAM_motif-chromVar_violinplots.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

###########################################
# now let's do something super jank
# each of the motifs we found have a set of associated peaks 
# let's just naively add all such peaks together, and that'll be a "motif score" for that cell

motif_peak_mapping <- all_data[["ATAC"]]@motifs@data
# keep only the relevant motifs
motif_peak_mapping <- motif_peak_mapping[,motifs]

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, c("dgCMatrix", "lgCMatrix")), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

activity_matrix <- sapply(1:ncol(motif_peak_mapping), function(i){
  peak_idx <- .nonzero_col(mat = motif_peak_mapping, col_idx = i, bool_value = F)
  Matrix::colSums(all_data[["ATAC"]]@data[peak_idx,])
})
colnames(activity_matrix) <- motif_names
activity_matrix <- Matrix::t(activity_matrix)

all_data[["motifActivity"]] <- Seurat::CreateAssayObject(activity_matrix)

Seurat::DefaultAssay(all_data) <- "motifActivity"
plot1 <- Seurat::VlnPlot(all_data, 
                         features = motif_names,
                         slot = "counts",
                         ncol = 6,
                         pt.size = 0.5)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6d/Writeup6d_DABTRAM_motifActivity_violinplots.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")


de_res2 <- Seurat::FindMarkers(
  object = all_data,
  test.use = "wilcox",
  ident.1 = paste0("day0_win_", treatment),
  ident.2 = paste0("day0_lose_", treatment),
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  verbose = F
)
de_res2
