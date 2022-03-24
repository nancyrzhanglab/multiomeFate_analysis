rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_processed2.RData")

conditions <- c("time10_cis", "time10_cocl2", "time10_dabtram")
main_vec <- c("CIS", "COCL2", "DABTRAM")
file_vec <- c("cis", "cocl2", "dabtram")
col_pal <- scales::hue_pal()(4)
base_col <- col_pal[1]; col_pal <- col_pal[-1]
### first, time0 and CIS

for(i in 1:length(conditions)){
  condition <- conditions[i]
  print(condition)
  
  all_data2 <- all_data
  keep_vec <- rep(0, ncol(all_data2))
  keep_vec[which(all_data2$original_dataset %in% c("time0", condition))] <- 1
  all_data2$keep <- keep_vec
  all_data2 <- subset(all_data2, keep == 1)
  
  Seurat::DefaultAssay(all_data2) <- "RNA"
  all_data2 <- Seurat::ScaleData(all_data2)
  all_data2 <- Seurat::RunPCA(all_data2, verbose = FALSE) 
  set.seed(10)
  all_data2 <- Seurat::RunUMAP(all_data2, dims = 1:20)
  
  Seurat::DefaultAssay(all_data2) <- "ATAC"
  all_data2 <- Signac::RunSVD(all_data2)
  set.seed(10)
  all_data2 <- Seurat::RunUMAP(object = all_data2, 
                               reduction = 'lsi', 
                               assay = 'ATAC', 
                               reduction.name = 'atac.umap', 
                               reduction.key = 'atacUMAP_',
                               dims = 2:25)
  
  col_vec <- c(base_col, col_pal[i])
  names(col_vec) <- c("time0", condition)
  plot1 <- Seurat::DimPlot(all_data2, reduction = "umap",
                           group.by = "original_dataset", label = TRUE,
                           repel = TRUE, label.size = 2.5,
                           cols = col_vec)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized, Time0 & Time10 ",
                                           main_vec[i], "),\n", length(all_data[["RNA"]]@var.features), 
                                           " genes, using 20 PCs"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_RNA_time0time10", 
                                    file_vec[i],".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  Seurat::DefaultAssay(all_data2) <- "RNA"
  plot1 <- Seurat::FeaturePlot(all_data2, 
                               reduction = "umap",
                               features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")),
                               ncol = 3)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_RNA_time0time10", 
                                    file_vec[i],"_jackpot1.png"),
                  plot1, device = "png", width = 12, height = 8, units = "in")
  
  plot1 <-Seurat::DimPlot(all_data2, reduction = "atac.umap",
                          group.by = "original_dataset", label = TRUE,
                          repel = TRUE, label.size = 2.5,
                          cols = col_vec)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (TF-IDF, Time0 & Time10 ",
                                           main_vec[i], "),\nusing 24 PCs"))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_ATAC_time0time10", 
                                    file_vec[i], ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
}


