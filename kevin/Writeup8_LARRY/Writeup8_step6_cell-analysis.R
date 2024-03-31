rm(list=ls())
library(Seurat)
library(multiomeFate)
library(biomaRt)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step3_fasttopics.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

seurat_object_safe <- seurat_object

treatment_vec <- as.character(sort(unique(seurat_object$time_celltype)))
day_early_vec <- treatment_vec[grep("^.*-4", treatment_vec)]
day_early <- "4"
treatment_vec <- treatment_vec[grep("^.*-6", treatment_vec)]

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
p1 <- p1 + ggplot2::ggtitle("S score\n(UMAP of RNA fasttopics)")
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_S-score_umap.png",
                p1, device = "png", width = 6.5, height = 5, units = "in")

p1 <- scCustomize::FeaturePlot_scCustom(seurat_object, 
                                        colors_use = list("red", "lightgray", "blue"),
                                        reduction = "umap", 
                                        na_cutoff = quantile(seurat_object$G2M.Score, probs = 0.05),
                                        na_color = "bisque",
                                        features = "G2M.Score")
p1 <- p1 + ggplot2::ggtitle("G2M score\n(UMAP of RNA fasttopics)")
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_G2M-score_umap.png",
                p1, device = "png", width = 6.5, height = 5, units = "in")

p1 <- Seurat::DimPlot(seurat_object, 
                      group.by = "Phase",
                      reduction = "umap")
p1 <- p1 + ggplot2::ggtitle("Cell-cycle phase\n(UMAP of RNA fasttopics)")
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_cellcycle-phase_umap.png",
                p1, device = "png", width = 6.5, height = 5, units = "in")

#############

cell_imputation_mat <- numeric(0)

for(treatment in treatment_vec){
  load(paste0("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_", treatment, "_from_day", day_early, "_postprocess.RData"))
  
  cell_imputation_mat <- cbind(cell_imputation_mat, cell_imputed_score)
}

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
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_percent_fate.png",
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
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_percent_fate_full.png",
                plot_all, device = "png", width = 14, height = 4, units = "in")


