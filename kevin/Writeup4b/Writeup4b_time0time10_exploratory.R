rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
source("seurat_helpers.R")

file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_02_arc_time0", "2022_02_arc_time10_CIS", 
                  "2022_02_arc_time10_COCL2", "2022_02_arc_time10_DABTRAM")

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

time0 <- create_seurat_object(file_folders[1])
time10_cis <- create_seurat_object(file_folders[2])
time10_cocl2 <- create_seurat_object(file_folders[3])
time10_dabtram <- create_seurat_object(file_folders[4])

starting_num_cells <- c(time0 = ncol(time0),
                        time10_cis = ncol(time10_cis),
                        time10_cocl2 = ncol(time10_cocl2),
                        time10_dabtram = ncol(time10_dabtram))

time0 <- normalize_seurat(time0)
time10_cis <- normalize_seurat(time10_cis)
time10_cocl2 <- normalize_seurat(time10_cocl2)
time10_dabtram <- normalize_seurat(time10_dabtram)

time0 <- qc_metrics(time0)
time10_cis <- qc_metrics(time10_cis)
time10_cocl2 <- qc_metrics(time10_cocl2)
time10_dabtram <- qc_metrics(time10_dabtram)

gene_intersection <- table(table(
  c(rownames(time0), rownames(time10_cis), 
    rownames(time10_cocl2), rownames(time10_dabtram))
))
peak_intersection <- table(table(
  c(rownames(time0[["ATAC"]]), rownames(time10_cis[["ATAC"]]), 
    rownames(time10_cocl2[["ATAC"]]), rownames(time10_dabtram[["ATAC"]]))
))

time0b <- time0; time0b[["ATAC"]] <- NULL
time10b_cis <- time10_cis; time10b_cis[["ATAC"]] <- NULL
time10b_cocl2 <- time10_cocl2; time10b_cocl2[["ATAC"]] <- NULL
time10b_dabtram <- time10_dabtram; time10b_dabtram[["ATAC"]] <- NULL

all_data <- merge(time0b, y = c(time10b_cis, time10b_cocl2, time10b_dabtram), 
                  add.cell.ids = c("time0", "time10_cis", "time10_cocl2", "time10_dabtram"), 
                  project = "All_Data", merge.data = T)
dataset_vec <- sapply(rownames(all_data@meta.data), function(x){
  tmp <- strsplit(x, split = "_")[[1]]
  if(length(tmp) == 2){
    return(tmp[1])
  } else {
    return(paste0(tmp[1:2], collapse = "_"))
  }
})
names(dataset_vec) <- NULL
all_data$original_dataset <- dataset_vec

Seurat::DefaultAssay(all_data) <- "RNA"
set.seed(10)
all_data <- Seurat::FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 500)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- intersect(jackpot_genes, rownames(all_data))
all_data[["RNA"]]@var.features <- unique(c(jackpot_genes, all_data[["RNA"]]@var.features))
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE) 
all_data <- Seurat::RunUMAP(all_data, dims = 1:25)

save(all_data, 
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_exploration_alldata_onlyGEX.RData")

#############################

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_jackpot1.png"),
                plot1, device = "png", width = 12, height = 8, units = "in")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(setdiff(jackpot_genes, c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E"))),
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_jackpot2.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")),
                             slot = "scale.data",
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_scaled_jackpot1.png"),
                plot1, device = "png", width = 12, height = 8, units = "in")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(setdiff(jackpot_genes, c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E"))),
                             slot = "scale.data",
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_umap_scaled_jackpot2.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

plot1 <-Seurat::DimPlot(all_data, reduction = "umap",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\n", length(all_data[["RNA"]]@var.features), " genes, using 30 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_500genes_umap_original_dataset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

############################

plot1 <- Seurat::VlnPlot(all_data, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
                         group.by = "original_dataset",
                         pt.size = 0.01,
                         ncol = 4)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_QC.png"),
                plot1, device = "png", width = 15, height = 6, units = "in")

############################

n <- ncol(all_data[["RNA"]])
tmp <- do.call(rbind, lapply(sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")), function(gene){
  cbind(all_data[["RNA"]]@counts[gene,], rep(gene, n), all_data$original_dataset)
}))
tmp <- as.data.frame(tmp)
colnames(tmp) <- c("value", "gene", "condition")
tmp$value <- as.numeric(tmp$value)
tmp$value[tmp$value > 10] <- 10

for(gene in c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")){
  tmp2 <- tmp[which(tmp$gene == gene),]
  print(max(tmp2$value))
  plot1 <- ggplot2::ggplot(tmp2, ggplot2::aes(y=condition, x=value,  fill=condition)) +
    ggridges::geom_density_ridges(alpha=0.6, stat="binline", bins=11,scale = 0.9) +
    ggridges::theme_ridges() +
    ggplot2::theme_gray() +
    Seurat::NoLegend() + 
    ggplot2::xlab(paste0("Observed count of ", gene)) +
    ggplot2::ylab("Condition")
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_jackpot_counts_", gene, ".png"),
                  plot1, device = "png", width = 10, height = 5, units = "in")
}


#############################

# try a umap with 2000 genes
all_data2 <- all_data
all_data2[["pca"]] <- NULL; all_data2[["umap"]] <- NULL
all_data2 <- Seurat::FindVariableFeatures(all_data2, selection.method = "vst", nfeatures = 2000)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- intersect(jackpot_genes, rownames(all_data2))
all_data2[["RNA"]]@var.features <- unique(c(jackpot_genes, all_data2[["RNA"]]@var.features))
all_data2 <- Seurat::ScaleData(all_data2)
all_data2 <- Seurat::RunPCA(all_data2, verbose = FALSE) 
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, dims = 1:30)

plot1 <-Seurat::DimPlot(all_data2, reduction = "umap",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\n", length(all_data2[["RNA"]]@var.features), " genes, using 30 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_2000genes_umap_original_dataset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#############################

# try a umap with 1000 genes
all_data2 <- all_data
all_data2[["pca"]] <- NULL; all_data2[["umap"]] <- NULL
all_data2 <- Seurat::FindVariableFeatures(all_data2, selection.method = "vst", nfeatures = 1000)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- intersect(jackpot_genes, rownames(all_data2))
all_data2[["RNA"]]@var.features <- unique(c(jackpot_genes, all_data2[["RNA"]]@var.features))
all_data2 <- Seurat::ScaleData(all_data2)
all_data2 <- Seurat::RunPCA(all_data2, verbose = FALSE) 
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, dims = 1:30)

plot1 <-Seurat::DimPlot(all_data2, reduction = "umap",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (LogNormalized),\n", length(all_data2[["RNA"]]@var.features), " genes, using 30 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_1000genes_umap_original_dataset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


