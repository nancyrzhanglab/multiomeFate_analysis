rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)
library(ordinal)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
metadata <- all_data@meta.data

treatment <- "DABTRAM"
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 5),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 25))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)

tier1_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier1_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier2_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier2_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier3_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier3_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
keep_vec <- rep(NA, ncol(all_data))
keep_vec[tier1_idx] <- paste0("3high_winner_", treatment)
keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
keep_vec[tier3_idx] <- paste0("1loser_", treatment)
names(keep_vec) <- colnames(all_data)

##############

load("../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM.RData")

tier_vec <- keep_vec[rownames(rna_mat)]
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "DABTRAM"
sd_vec <- sparseMatrixStats::colSds(rna_mat)
if(any(sd_vec <= 1e-6)){
  rna_mat <- rna_mat[,sd_vec >= 1e-6,drop = F]
}
gene_vec <- sort(intersect(names(chr_peak_list), colnames(rna_mat)))

n <- nrow(rna_mat)
y <- multiomeFate:::form_onehot_classification_mat(tier_vec)

spca_res_list <- vector("list", length = length(gene_vec))
names(spca_res_list) <- gene_vec
for(i in 1:length(gene_vec)){
  gene <- gene_vec[i]
  print(paste0(gene, ": ", i, " out of ", length(gene_vec)))
  chr_peak_mat <- chr_peak_list[[gene]]
  if(!is.matrix(chr_peak_mat)) {
    chr_peak_mat <- matrix(chr_peak_mat, nrow = length(chr_peak_mat), ncol = 1)
    colnames(chr_peak_mat) <- paste0(gene, ":ATAC")
  }
  tmp <- cbind(rna_mat[,gene], chr_peak_mat)
  colnames(tmp)[1] <- paste0(gene, ":RNA")
  tmp <- scale(tmp)
  spca_res_list[[gene]] <- multiomeFate:::supervised_pca(x = tmp, y = y)
}

save(spca_res_list, date_of_run, session_info, tab_mat,
     metadata, tier_vec, rna_mat,
     chr_peak_list,
     file = "../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_spca.RData")

# head(spca_res_list[["FN1"]]$U)
# head(spca_res_list[["FN1"]]$dimred)
# 
# percent_mat <- sapply(spca_res_list, function(zz){
#   tmp <- zz$U
#   (tmp[1,]^2)*100
# })
# round(percent_mat)

##############################

source("ordinal_functions.R")

y_vec <- rep(NA, length(tier_vec))
idx_list <- vector("list", length = 3)
idx_list[[1]] <- which(tier_vec == paste0("1loser_", treatment))
idx_list[[2]] <- which(tier_vec == paste0("2mid_winner_", treatment))
idx_list[[3]] <- which(tier_vec == paste0("3high_winner_", treatment))
for(i in 1:3){
  y_vec[idx_list[[i]]] <- i
}
y_vec <- as.factor(y_vec)

cv_score_vec <- rep(NA, length(spca_res_list))
names(cv_score_vec) <- names(spca_res_list)
for(i in 1:length(spca_res_list)){
  gene <- names(spca_res_list)[i]
  print(paste0(gene, ": ", i, " out of ", length(gene_vec)))
  set.seed(10)
  
  x_mat <- Re(spca_res_list[[gene]]$dimred)
  cv_score_vec[gene] <- .five_fold_cv(x_mat, y_vec)
}
table(is.na(cv_score_vec))
round(quantile(100*cv_score_vec, na.rm = T))

for(i in 1:length(spca_res_list)){
  spca_res_list[[i]]$U <- Re(spca_res_list[[i]]$U)
  spca_res_list[[i]]$dimred <- Re(spca_res_list[[i]]$dimred)
}


save(spca_res_list, date_of_run, session_info, tab_mat,
     metadata, tier_vec, rna_mat,
     chr_peak_list, cv_score_vec,
     file = "../../../../out/kevin/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_spca.RData")

########

# null_dist_list <- vector("list", length = length(spca_res_list))
# names(null_dist_list) <- names(spca_res_list)
# for(gene in names(spca_res_list)){
#   print(gene)
#   x_mat <- spca_res_list[[gene]]$dimred
#   
#   set.seed(10)
#   null_dist_list[[gene]] <- .permutation_null_score(x_mat, y_vec, trials = 100)
# }
# 
# p_val <- sapply(1:length(null_dist_list), function(i){
#   length(which(cv_score_vec[i] <= null_dist_list[[i]]))/length(null_dist_list[[i]])
# })

##########

percent_rna <- sapply(gene_vec, function(gene){
  100*Re(spca_res_list[[gene]]$U[1,1])^2
})


source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec2 <- sort(unique(c(unlist(keygenes), keygenes_csc)))

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

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_keygenes-and-chrpeak_DABTRAM_rnapercentage_cvscore2.png"),
                p1, device = "png", width = 10, height = 10, units = "in")

################################

# gene_vec <- intersect(names(cv_score_vec)[which(cv_score_vec >= 0.425)],
#                       names(percent_rna)[which(percent_rna <= 0.4)]) 
gene_vec <- names(cv_score_vec)[which(cv_score_vec >= 0.5)]
gene_vec <- sort(gene_vec)

tmp_list <- lapply(gene_vec, function(gene){
  Re(spca_res_list[[gene]]$dimred)
})
spca_mat <- do.call(cbind, tmp_list)

set.seed(10)
seurat_obj <- Seurat::CreateSeuratObject(counts = t(spca_mat))
seurat_obj[["RNA"]]@data <- seurat_obj[["RNA"]]@counts
seurat_obj[["RNA"]]@scale.data <- as.matrix(seurat_obj[["RNA"]]@counts)
seurat_obj[["RNA"]]@var.features <- rownames(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = F)

zz <- seurat_obj[["pca"]]@cell.embeddings
sing_vec <- sapply(1:ncol(zz), function(j){sqrt(sum(zz[,j]^2))})
round(100*diff(sing_vec)/sing_vec[-1])
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:50)
seurat_obj$tier_vec <- tier_vec

col_palette <- c("gray", "blue", "red")
names(col_palette) <- sort(unique(tier_vec))
p1 <- Seurat::DimPlot(seurat_obj, reduction = "umap",
                      cols = col_palette,
                      group.by = "tier_vec")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_spca_day10.png"),
                p1, device = "png", width = 7, height = 5, units = "in")



