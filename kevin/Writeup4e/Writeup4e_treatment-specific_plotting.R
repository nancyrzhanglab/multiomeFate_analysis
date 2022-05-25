rm(list=ls())

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["Lineage"]] <- NULL
all_data[["spliced"]] <- NULL
all_data[["unspliced"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["pca"]] <- NULL
all_data[["umap"]] <- NULL
all_data[["lsi"]] <- NULL
all_data[["atac.umap"]] <- NULL
all_data[["activityPCA"]] <- NULL
all_data[["activity.umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["saverumap"]] <- NULL

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
for(treatment in treatment_vec){
  print(treatment)
  
  keep_vec <- rep(0, ncol(all_data))
  keep_vec[c(which(all_data$dataset == "day0"), grep(treatment, all_data$dataset))] <- 1
  all_data$keep <- keep_vec
  all_data_subset <- subset(all_data, keep == 1)
  
  # ATAC
  Seurat::DefaultAssay(all_data_subset) <- "ATAC"
  all_data_subset <- Signac::RunTFIDF(all_data_subset)
  all_data_subset <- Signac::RunSVD(all_data_subset)
  set.seed(10)
  all_data_subset <- Seurat::RunUMAP(object = all_data_subset, 
                                     reduction = 'lsi', dims = 2:50,
                                     reduction.name = 'atac.umap')
  
  # gene activity
  Seurat::DefaultAssay(all_data_subset) <- "geneActivity"
  all_data_subset <- Seurat::RunPCA(all_data_subset, verbose = FALSE,
                                    reduction.name = "activityPCA") 
  set.seed(10)
  all_data_subset <- Seurat::RunUMAP(all_data_subset, dims = 2:50,
                                     reduction = "activityPCA", 
                                     reduction.name = "activity.umap")
  
  
  # SAVER
  Seurat::DefaultAssay(all_data_subset) <- "Saver"
  all_data_subset <- Seurat::RunPCA(all_data_subset, verbose = FALSE,
                                    reduction.name = "saverpca") 
  set.seed(10)
  all_data_subset <- Seurat::RunUMAP(all_data_subset, dims = 1:50,
                                     reduction = "saverpca",
                                     reduction.name = "saverumap")
  
  # plotting
  Seurat::DefaultAssay(all_data_subset) <- "ATAC"
  plot1 <- Seurat::DimPlot(all_data_subset, 
                           reduction = "atac.umap", 
                           group.by = "dataset", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (TF-IDF): ", treatment, ",\nusing 49 PCs"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_atac_umap-", treatment, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  Seurat::DefaultAssay(all_data_subset) <- "geneActivity"
  plot1 <- Seurat::DimPlot(all_data_subset, 
                           reduction = "activity.umap", 
                           group.by = "dataset", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (Gene activity): ", treatment, ",\nusing 49 PCs"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_atac-activity_umap-", treatment, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  Seurat::DefaultAssay(all_data_subset) <- "Saver"
  plot1 <-Seurat::DimPlot(all_data_subset, reduction = "saverumap",
                          group.by = "dataset", label = TRUE,
                          repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (Saver) : ", treatment,"\n,", length(all_data_subset[["Saver"]]@var.features), " genes, using 50 PCs"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4e/Writeup4e_saver_umap-", treatment, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}