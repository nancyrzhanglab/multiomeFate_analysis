rm(list=ls())
load("../../../../out/kevin/Writeup4c/Writeup4c_timeAll_peakmerging_simplified.RData")
library(Seurat)
library(Signac)

treatment_vec <- c("CIS", "COCL2", "DABTRAM", "day10", "week5")
keep_list <- lapply(treatment_vec, function(treatment){
  keep_vec <- rep(0, ncol(all_data))
  if(treatment %in% c("CIS", "COCL2", "DABTRAM")){
    keep_vec[which(all_data$original_dataset %in% c("day0", paste0("day10_",treatment), paste0("week5_",treatment)))] <- 1
  } else {
    keep_vec[c(grep("day0", all_data$original_dataset), grep(treatment, all_data$original_dataset))] <- 1
  }
  keep_vec
})
names(keep_list) <- treatment_vec

for(i in 4:5){
  seurat_obj <- all_data
  treatment <- names(keep_list)[i]
  seurat_obj$keep <- keep_list[[i]]
  seurat_obj <- subset(seurat_obj, keep == 1)
  seurat_obj$pca <- NULL
  seurat_obj$lsi <- NULL
  seurat_obj$umap <- NULL
  seurat_obj$adt.umap <- NULL
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = FALSE) 
  set.seed(10)
  rna_npcs <- 30
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:rna_npcs)
  
  Seurat::DefaultAssay(seurat_obj) <- "ATAC"
  set.seed(10)
  seurat_obj <- Signac::RunSVD(seurat_obj)
  set.seed(10)
  atac_npcs <- 35
  seurat_obj <- Seurat::RunUMAP(object = seurat_obj, 
                              reduction = 'lsi', dims = 2:atac_npcs,
                              reduction.name = 'adt.umap')
  
  ##########
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  plot1 <-Seurat::DimPlot(seurat_obj, reduction = "umap",
                          group.by = "dataset", label = TRUE,
                          repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\n", length(all_data[["RNA"]]@var.features), " genes, using ", rna_npcs, " PCs: ", treatment))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4c/Writeup4c_rna_umap_", treatment, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  
  Seurat::DefaultAssay(seurat_obj) <- "ATAC"
  plot1 <- Seurat::DimPlot(seurat_obj, 
                           reduction = "adt.umap", 
                           group.by = "dataset", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (TF-IDF),\nusing ", atac_npcs-1, " PCs: ", treatment))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4c/Writeup4c_atac_umap_", treatment, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  plot1 <- Seurat::FeaturePlot(seurat_obj, 
                               reduction = "umap",
                               slot = "scale.data",
                               features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")),
                               ncol = 3)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4c/Writeup4c_rna_umap_jackpot1_", treatment, ".png"),
                  plot1, device = "png", width = 12, height = 8, units = "in")
  
}

##################


