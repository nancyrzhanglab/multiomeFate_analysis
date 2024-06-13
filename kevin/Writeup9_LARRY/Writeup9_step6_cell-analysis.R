rm(list=ls())
library(patchwork)
library(ggplot2)
library(ggtern)
library(Seurat)
library(multiomeFate)
library(biomaRt)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step3_fasttopics.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

seurat_object_safe <- seurat_object

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
day_early <- "4"
day_later <- "6"
day_early_vec <- treatment_vec[grep(paste0("^.*-", day_early), treatment_vec)]
treatment_vec <- treatment_vec[grep(paste0("^.*-", day_later), treatment_vec)]

# https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# https://github.com/satijalab/seurat/issues/2493
# https://support.bioconductor.org/p/129636/ and https://support.bioconductor.org/p/9144001/
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol"), 
                            filters = "hgnc_symbol", 
                            values = x , 
                            mart = human, 
                            attributesL = c("mgi_symbol"), 
                            martL = mouse, 
                            uniqueRows = T)
  
  humanx <- unique(genesV2[, 2])
  
  return(humanx)
}

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
m.s.genes <- convertHumanGeneList(s.genes)
m.g2m.genes <- convertHumanGeneList(g2m.genes)
seurat_object <- Seurat::CellCycleScoring(seurat_object, 
                                          s.features = m.s.genes, 
                                          g2m.features = m.g2m.genes)
seurat_object_safe <- seurat_object

#############

seurat_object <- seurat_object_safe
keep_vec <- rep(FALSE, ncol(seurat_object))
idx <- which(seurat_object$time_celltype %in% day_early_vec)
keep_vec[idx] <- TRUE
seurat_object$keep <- keep_vec
seurat_object <- subset(seurat_object, keep == TRUE)

p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        reduction = "umap", 
                                        na_cutoff = quantile(seurat_object$S.Score, probs = 0.05),
                                        na_color = "bisque",
                                        features = "S.Score")
p1 <- p1 + ggplot2::ggtitle("S score\n(UMAP of RNA)")
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_S-score_", day_early, "_umap.png"),
                p1, device = "png", width = 6.5, height = 5, units = "in")
# https://github.com/tidyverse/ggplot2/issues/5612 bizarre error of "Error in Ops.data.frame(guide_loc, panel_loc) :" fixed using this

p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        reduction = "umap", 
                                        na_cutoff = quantile(seurat_object$G2M.Score, probs = 0.05),
                                        na_color = "bisque",
                                        features = "G2M.Score")
p1 <- p1 + ggplot2::ggtitle("G2M score\n(UMAP of RNA)")
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_G2M-score_", day_early, "_umap.png"),
                p1, device = "png", width = 6.5, height = 5, units = "in")

p1 <- Seurat::DimPlot(seurat_object, 
                      group.by = "Phase",
                      reduction = "umap")
p1 <- p1 + ggplot2::ggtitle("Cell-cycle phase\n(UMAP of RNA)")
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_cellcycle-phase_", day_early, "_umap.png"),
                p1, device = "png", width = 6.5, height = 5, units = "in")

#############

cell_imputation_mat <- numeric(0)

for(treatment in treatment_vec){
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment, "_from_day", day_early, "_postprocess.RData"))
  
  cell_imputation_mat <- cbind(cell_imputation_mat, cell_imputed_score)
}
colnames(cell_imputation_mat) <- treatment_vec

cell_imputation_mat2 <- 10^cell_imputation_mat
# cell_imputation_mat2 <- pmax(cell_imputation_mat2, 10^(-1.5))
n <- nrow(cell_imputation_mat2)
for(i in 1:n){
  tmp <- cell_imputation_mat2[i,]
  if(sum(tmp) <= 0.01){
    cell_imputation_mat2[i,] <- NA
  } else {
    cell_imputation_mat2[i,] <- tmp/sum(tmp)
  }
}
colnames(cell_imputation_mat2) <- paste0("percent_fate_", treatment_vec)
idx <- unique(unlist(apply(cell_imputation_mat2, 2, function(x){which(is.na(x))})))
if(length(idx) > 0) cell_imputation_mat2 <- cell_imputation_mat2[-idx,,drop = FALSE]


for(j in 1:3){
  vec <- rep(NA, length = length(SeuratObject::Cells(seurat_object)))
  names(vec) <- SeuratObject::Cells(seurat_object)
  
  vec[rownames(cell_imputation_mat2)] <- cell_imputation_mat2[,j]
  seurat_object@meta.data[,colnames(cell_imputation_mat2)[j]] <- vec
}

plot_list <- lapply(treatment_vec, function(treatment){
  colname_val <- paste0("percent_fate_", treatment)
  vec <- seurat_object@meta.data[,colname_val]
  
  p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                          colors_use = list("red", "lightgray", "blue"),
                                          reduction = "umap", 
                                          na_cutoff = quantile(vec, probs = 0.05, na.rm = TRUE),
                                          na_color = "bisque",
                                          features = colname_val)
  p1 <- p1 + ggplot2::ggtitle(paste0("Percent fate to ", treatment, "\n(UMAP of RNA fasttopics)"))
  p1
})

plot_all <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_percent_fate_", day_early, ".png"),
                plot_all, device = "png", width = 14, height = 4, units = "in")

# without thresholding
cell_imputation_mat2 <- 10^cell_imputation_mat
n <- nrow(cell_imputation_mat2)
for(i in 1:n){
  tmp <- cell_imputation_mat2[i,]
  cell_imputation_mat2[i,] <- tmp/sum(tmp)
}
colnames(cell_imputation_mat2) <- paste0("percent_fate_", treatment_vec)

for(j in 1:3){
  vec <- rep(NA, length = length(SeuratObject::Cells(seurat_object)))
  names(vec) <- SeuratObject::Cells(seurat_object)
  
  vec[rownames(cell_imputation_mat2)] <- cell_imputation_mat2[,j]
  seurat_object@meta.data[,colnames(cell_imputation_mat2)[j]] <- vec
}

plot_list <- lapply(treatment_vec, function(treatment){
  colname_val <- paste0("percent_fate_", treatment)
  vec <- seurat_object@meta.data[,colname_val]
  
  p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                          colors_use = list("red", "lightgray", "blue"),
                                          reduction = "umap", 
                                          na_cutoff = quantile(vec, probs = 0.05, na.rm = TRUE),
                                          na_color = "bisque",
                                          features = colname_val)
  p1 <- p1 + ggplot2::ggtitle(paste0("Percent fate to ", treatment, "\n(UMAP of RNA fasttopics)"))
  p1
})

plot_all <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_percent_fate_full_", day_early, ".png"),
                plot_all, device = "png", width = 14, height = 4, units = "in")

#############################

# pairs plot of the genes
rna_mat <- SeuratObject::LayerData(seurat_object,
                                   data = "data",
                                   assay = "RNA",
                                   features = Seurat::VariableFeatures(seurat_object))
gene_corr_list <- vector("list", length = length(treatment_vec))
names(gene_corr_list) <- treatment_vec

for(treatment in treatment_vec){
  rna_mat2 <- Matrix::t(rna_mat[,rownames(cell_imputation_mat)])
  corr_vec <- sapply(1:ncol(rna_mat2), function(j){
    stats::cor(cell_imputation_mat[,treatment], rna_mat2[,j])
  })
  names(corr_vec) <- colnames(rna_mat2)
  corr_vec[is.na(corr_vec)] <- 0
  gene_corr_list[[treatment]] <- corr_vec
}

df <- data.frame(do.call(cbind, gene_corr_list))
p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_growth-potential_", day_early, "_pairs_gene-corr.png"),
                p1, device = "png", width = 8, height = 8, units = "in")

####

# from supplementary table 2 in https://www.science.org/doi/full/10.1126/science.aaw3381
marker_list <- list(
  Monocyte = c("Ms4a6d", "Fabp5", "Ctss", "Ms4a6c", "Tgfbi", "Olmf1", "Csfl1r", "Ccr2", "Klf4", "F13a1"),
  Neutrophil = c("S100a9", "Itgb21", "Elane", "Fcnb", "Mpo", "Prtn3", "S100a6", "S100a8", "Lcn2", "Lrg1") 
)

for(i in 1:length(marker_list)){
  celltype <- names(marker_list)[i]
  print(celltype)
  
  marker_vec <- marker_list[[i]]
  marker_vec <- intersect(marker_vec, names(gene_corr_list[[1]]))
  print(paste0("Number of markers: ", length(marker_vec)))
  
  df <- data.frame(do.call(cbind, gene_corr_list))
  col_vec <- rep("coral1", nrow(df)); names(col_vec) <- rownames(df)
  cex_vec <- rep(1, nrow(df)); names(cex_vec) <- rownames(df)
  
  gene_idx <- which(rownames(df) %in% marker_vec)
  col_vec[gene_idx] <- "black"
  cex_vec[gene_idx] <- 2
  gene_ordering <- c(setdiff(1:nrow(df), gene_idx), gene_idx)
  
  df <- df[gene_ordering,]
  col_vec <- col_vec[gene_ordering]
  cex_vec <- cex_vec[gene_ordering]
  
  png(paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_growth-potential_", day_early, "_pairs_gene-corr_highlight_",
             celltype,
             ".png"),
      height = 2500, width = 2500, units = "px", res = 300)
  graphics::pairs(df, col = col_vec, pch = 16, cex = cex_vec,
                  main = celltype)
  graphics.off()
}

#############################

# pairs plot of the genes, now based on their coefficients

df <- numeric(0)
load(paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment_vec[1], "_from_day", day_early, "_postprocess.RData"))
vec <- fit$coefficient_vec
gene_names <- names(vec)
gene_names <- gene_names[!gene_names %in% "Intercept"]

for(treatment in treatment_vec){
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment, "_from_day", day_early, "_postprocess.RData"))
  
  vec <- fit$coefficient_vec
  vec <- vec[gene_names]
  
  df <- cbind(df, vec)
}
colnames(df) <- treatment_vec

df <- data.frame(df)

p1 <- GGally::ggpairs(df, 
                      mapping = ggplot2::aes(color = "coral1"),
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_growth-potential_", day_early, "_pairs_gene-corr_coefficients.png"),
                p1, device = "png", width = 8, height = 8, units = "in")

####

# from supplementary table 2 in https://www.science.org/doi/full/10.1126/science.aaw3381
marker_list <- list(
  Monocyte = c("Ms4a6d", "Fabp5", "Ctss", "Ms4a6c", "Tgfbi", "Olmf1", "Csfl1r", "Ccr2", "Klf4", "F13a1"),
  Neutrophil = c("S100a9", "Itgb21", "Elane", "Fcnb", "Mpo", "Prtn3", "S100a6", "S100a8", "Lcn2", "Lrg1") 
)

for(i in 1:length(marker_list)){
  celltype <- names(marker_list)[i]
  print(celltype)
  
  marker_vec <- marker_list[[i]]
  marker_vec <- intersect(marker_vec, gene_names)
  print(paste0("Number of markers: ", length(marker_vec)))
  
  df2 <- df
  col_vec <- rep("coral1", nrow(df2)); names(col_vec) <- rownames(df2)
  cex_vec <- rep(1, nrow(df2)); names(cex_vec) <- rownames(df2)
  
  gene_idx <- which(rownames(df2) %in% marker_vec)
  col_vec[gene_idx] <- "black"
  cex_vec[gene_idx] <- 2
  gene_ordering <- c(setdiff(1:nrow(df2), gene_idx), gene_idx)
  
  df2 <- df2[gene_ordering,]
  col_vec <- col_vec[gene_ordering]
  cex_vec <- cex_vec[gene_ordering]
  
  png(paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_growth-potential_", day_early, "_pairs_gene-corr_highlight_",
             celltype,
             "_coefficients.png"),
      height = 2500, width = 2500, units = "px", res = 300)
  graphics::pairs(df2, col = col_vec, pch = 16, cex = cex_vec,
                  main = celltype)
  graphics.off()
}

#################################

# pairs plot of the cells
cell_imputation_mat <- numeric(0)

for(treatment in treatment_vec){
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup9/Writeup9_", treatment, "_from_day", day_early, "_postprocess.RData"))
  
  cell_imputation_mat <- cbind(cell_imputation_mat, cell_imputed_score)
}
colnames(cell_imputation_mat) <- treatment_vec

df <- data.frame(cell_imputation_mat)
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_growth-potential_", day_early, "_pairs_cell-corr_logscale.png"),
                p1, device = "png", width = 8, height = 8, units = "in")


df <- data.frame(10^cell_imputation_mat)
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup9/Writeup9_growth-potential_", day_early, "_pairs_cell-corr.png"),
                p1, device = "png", width = 8, height = 8, units = "in")


