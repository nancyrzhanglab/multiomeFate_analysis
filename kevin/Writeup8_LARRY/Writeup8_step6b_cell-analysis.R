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
day_early_vec <- treatment_vec[grep("^.*-2", treatment_vec)]
day_early <- "2"
treatment_vec <- treatment_vec[grep("^.*-4", treatment_vec)]

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
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_S-score_", day_early, "_umap.png"),
                p1, device = "png", width = 6.5, height = 5, units = "in")
# https://github.com/tidyverse/ggplot2/issues/5612 bizarre error of "Error in Ops.data.frame(guide_loc, panel_loc) :" fixed using this

p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        reduction = "umap", 
                                        na_cutoff = quantile(seurat_object$G2M.Score, probs = 0.05),
                                        na_color = "bisque",
                                        features = "G2M.Score")
p1 <- p1 + ggplot2::ggtitle("G2M score\n(UMAP of RNA)")
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_G2M-score_", day_early, "_umap.png"),
                p1, device = "png", width = 6.5, height = 5, units = "in")

p1 <- Seurat::DimPlot(seurat_object, 
                      group.by = "Phase",
                      reduction = "umap")
p1 <- p1 + ggplot2::ggtitle("Cell-cycle phase\n(UMAP of RNA)")
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_cellcycle-phase_", day_early, "_umap.png"),
                p1, device = "png", width = 6.5, height = 5, units = "in")

#############

cell_imputation_mat <- numeric(0)

for(treatment in treatment_vec){
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_", treatment, "_from_day", day_early, "_postprocess.RData"))
  
  cell_imputation_mat <- cbind(cell_imputation_mat, cell_imputed_score)
}
colnames(cell_imputation_mat) <- treatment_vec

cell_imputation_mat2 <- 10^cell_imputation_mat
cell_imputation_mat2 <- pmax(cell_imputation_mat2, 10^(-1.5))
n <- nrow(cell_imputation_mat2)
for(i in 1:n){
  tmp <- cell_imputation_mat2[i,]
  if(sum(tmp) <= 0.1){
    cell_imputation_mat2[i,] <- NA
  } else {
    cell_imputation_mat2[i,] <- tmp/sum(tmp)
  }
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
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_percent_fate_", day_early, ".png"),
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
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_percent_fate_full_", day_early, ".png"),
                plot_all, device = "png", width = 14, height = 4, units = "in")

######################################

# pairs plot of the cells

cell_imputation_mat2 <- 10^cell_imputation_mat
cell_imputation_mat2 <- pmax(cell_imputation_mat2, 10^(-1.5))
n <- nrow(cell_imputation_mat2)
for(i in 1:n){
  tmp <- cell_imputation_mat2[i,]
  if(sum(tmp) <= 0.1){
    cell_imputation_mat2[i,] <- NA
  } else {
    cell_imputation_mat2[i,] <- tmp/sum(tmp)
  }
}
colnames(cell_imputation_mat2) <- paste0("percent_fate_", treatment_vec)
celltype_vec <- as.character(seurat_object$time_celltype[rownames(cell_imputation_mat2)])
celltype_vec <- sapply(celltype_vec, function(val){
  strsplit(val, split = "-")[[1]][1]
})
names(celltype_vec) <- NULL

# https://cran.r-project.org/web/packages/Ternary/vignettes/Ternary.html ?
# https://www.marvinschmitt.com/blog/ggsimplex-prerelease/
# https://cran.r-project.org/web/packages/ggtern/index.html
# https://rpubs.com/KDVdecisions/triadtutorial1
# http://www.ggtern.com/d/2.2.2/ggsave.html

# add jitter
cell_imputation_mat3 <- cell_imputation_mat2
set.seed(10)
n <- nrow(cell_imputation_mat3)
for(i in 1:n){
  if(any(is.na(cell_imputation_mat3[i,]))) next()
  cell_imputation_mat3[i,] <- cell_imputation_mat3[i,] + stats::runif(3, min = 0, max = 0.1)
  cell_imputation_mat3[i,] <- cell_imputation_mat3[i,]/sum(cell_imputation_mat3[i,])
}

cell_imputation_mat3 <- as.data.frame(cell_imputation_mat3)
colnames(cell_imputation_mat3) <- c("Monocyte", "Neutrophil", "Undifferentiated")
cell_imputation_mat3 <- cbind(cell_imputation_mat3, celltype_vec)
colnames(cell_imputation_mat3)[4] <- "celltype"
cell_imputation_mat3$celltype <- factor(cell_imputation_mat3$celltype)
idx <- which(is.na(cell_imputation_mat3[,1]))
cell_imputation_mat3 <- cell_imputation_mat3[-idx,]

color_palette <- c("blue3", "coral2", "gray50")
names(color_palette) <- paste0(c("Monocyte", "Neutrophil", "Undifferentiated"))
tab_vec <- table(cell_imputation_mat3$celltype)
tab_vec <- sort(tab_vec, decreasing = TRUE)
row_idx <- unlist(lapply(names(tab_vec), function(val){
  which(cell_imputation_mat3$celltype == val)
}))
cell_imputation_mat3 <- cell_imputation_mat3[row_idx,]

plot1 <- ggtern::ggtern(data = cell_imputation_mat3,
                        mapping = ggplot2::aes(x = Monocyte, 
                                               y = Neutrophil, 
                                               z = Undifferentiated,
                                               color = celltype)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = color_palette) +
  ggtern::theme_showarrows() + 
  ggplot2::labs(x = "Monocyte", 
                y = "Neutrophil",
                z = "Undifferentiated",
                title = paste0("Percent fate for day ", day_early))
ggtern::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_percent-fate_", day_early, "_cell-pairs.png"),
               plot1, 
               device = "png", 
               width = 10, 
               height = 5, 
               units = "in")

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
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_growth-potential_", day_early, "_pairs_gene-corr.png"),
                p1, device = "png", width = 8, height = 8, units = "in")


