rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")
load("../../../../out/kevin/Writeup6g/gene_vec.RData")

tmp <- sapply(granges_list, length)
granges_list <- granges_list[which(tmp > 0)]
granges_list <- granges_list[order(names(granges_list))]

# extend each one
for(i in 1:length(granges_list)){
  if(i %% floor(length(granges_list)/10) == 0) cat('*')
  # from https://github.com/stuart-lab/signac/blob/master/R/utilities.R#L1271
  granges_list[[i]] <- suppressWarnings(Signac::Extend(granges_list[[i]],
                                      upstream = 5000,
                                      downstream = 5000))
}

all_data[["ATAC"]]@motifs@data[1:5,1:5] # peak-by-motif for where peaks a motif is located within
length(all_data[["ATAC"]]@motifs@motif.names)
all_data[["ATAC"]]@motifs@positions[[1]] # the specific location a motif is detected to be in (observe: the range is quite short)
length(all_data[["ATAC"]]@motifs@positions[[1]])
length(multiomeFate:::.nonzero_col(all_data[["ATAC"]]@motifs@data, col_idx = 1, bool_value = F)) # smaller since multiple motif-binding might occur within a peak

all_data[["chromvar"]]@data[1:5,1:5] # motif-by-cell for the chromVar score of each cell

# let's just focus on the "important" genes for now
load("../../../../out/kevin/Writeup6g/custom_peakAggregation_tmp.RData")
len_vec <- sapply(result_list, length)
result_list <- result_list[which(len_vec > 0)]

granges_list <- granges_list[names(result_list)]
# for each gene, figure out what the peak-indices are
peakname_vec <- rownames(all_data[["ATAC"]]@motifs@data)
peak_idx_list <- vector("list", length(result_list))
names(peak_idx_list) <- names(result_list)
for(i in 1:length(result_list)){
  if(i %% floor(length(result_list)/10) == 0) cat('*')
  
  mat <- result_list[[i]]
  colname_vec <- colnames(mat)
  tmp <- strsplit(colname_vec, split = ":")
  gene_peakname_vec <- sapply(1:length(tmp), function(z){
    paste0(tmp[[z]][-1], collapse = "-")
  })
  peak_idx_list[[i]] <- which(peakname_vec %in% gene_peakname_vec)
}
quantile(sapply(peak_idx_list, length)) # some might be "smaller than expected" (or 0), since when I processed the peaks, I truncated the ends that were on the boundary of the +/- 5000bp

# now let's get that vector of important motifs
treatment <- "DABTRAM"

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
length(tier1_lineages); length(tier3_lineages)

tier3_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier3_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier1_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier1_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
ident_vec <- rep(NA, ncol(all_data))
names(ident_vec) <- colnames(all_data)
ident_vec[tier3_idx] <- "winner"
ident_vec[tier1_idx] <- "loser"
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"

set.seed(10)
Seurat::DefaultAssay(all_data) <- "chromvar"
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = "winner",
  ident.2 = "loser",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  verbose = F
)

idx <- which(de_res[,"p_val_adj"] <= 1e-2)
if(length(idx) == 0) idx <- 1
print(length(idx))
de_res <- de_res[idx,]
motifs <- rownames(de_res)
motif_names <- Signac::GetMotifData(object = all_data,
                                 assay = "ATAC",
                                 slot = "motif.names")[motifs]
cbind(motif_names, de_res[,"avg_diff"])

binding_mat <- matrix(NA, nrow = length(peak_idx_list), ncol = length(motifs))
rownames(binding_mat) <- names(peak_idx_list)
colnames(binding_mat) <- motifs

for(i in 1:length(motifs)){
  print(paste0(i, " out of ", length(motifs)))
  
  motif <- motifs[i]
  motif_idx <- which(colnames(all_data[["ATAC"]]@motifs@data) == motif)
  peak_idx <- multiomeFate:::.nonzero_col(all_data[["ATAC"]]@motifs@data, col_idx = motif_idx, bool_value = F)
  
  # run through peak_idx_list
  binding_mat[,i] <- sapply(peak_idx_list, function(vec){
    length(which(vec %in% peak_idx))
  })
}
quantile(binding_mat)
length(which(binding_mat>10))
quantile(apply(binding_mat,2,function(x){length(which(x>20))}))

num_peak_per_gene <- sapply(peak_idx_list, length)
num_peak_per_gene <- pmax(num_peak_per_gene, 1)
binding_mat2 <- diag(1/num_peak_per_gene) %*% binding_mat 

##############

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec2 <- sort(unique(c(unlist(keygenes), keygenes_csc)))
gene_vec <- rownames(binding_mat)

pdf("../../../../out/figures/Writeup6g/Writeup6g_chromVar_motif-locations_DABTRAM.pdf", 
    onefile = T, width = 5, height = 5)
for(j in 1:ncol(binding_mat)){
  print(j)
  
  df <- data.frame(number_peaks = binding_mat[,j],
                   percent_peaks = binding_mat2[,j],
                   gene = gene_vec)
  label_vec <- rep(F, length(gene_vec))
  names(label_vec) <- gene_vec
  idx <- which(gene_vec %in% gene_vec2)
  idx2 <- unique(c(which(binding_mat[,j] >= max(quantile(binding_mat[,j], prob = 0.9),1)),
                   which(binding_mat2[,j] >= max(quantile(binding_mat2[,j], prob = 0.9),0.1))))
  label_vec[intersect(idx,idx2)] <- T
  df$labeling = label_vec
  df <- df[c(which(!df[,"labeling"]), which(df[,"labeling"])),]
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = number_peaks, y = percent_peaks))
  p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
  p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
  p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == TRUE),
                                      ggplot2::aes(label = gene, color = labeling),
                                      size = 2,
                                      max.overlaps = 50)
  p1 <- p1 + ggplot2::ggtitle(paste0(motifs[j], ": ", motif_names[j], 
                                     "\n-Log10 pval=", round(-log10(de_res[motifs[j],"p_val"]), 2),
                                     ", value=", round(de_res[motifs[j],"avg_diff"], 2))) + Seurat::NoLegend()
  print(p1)
}

dev.off()