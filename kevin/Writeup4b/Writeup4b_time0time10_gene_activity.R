rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
source("seurat_helpers.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

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

###################3

# https://satijalab.org/signac/reference/geneactivity
set.seed(10)
Seurat::DefaultAssay(time0) <- "ATAC"
time0_activity <- Signac::GeneActivity(time0, 
                                       extend.downstream = 5000,
                                       extend.upstream = 5000,
                                       verbose = T)
set.seed(10)
Seurat::DefaultAssay(time10_cis) <- "ATAC"
time10_cis_activity <- Signac::GeneActivity(time10_cis, 
                                       extend.downstream = 5000,
                                       extend.upstream = 5000,
                                       verbose = T)
set.seed(10)
Seurat::DefaultAssay(time10_cocl2) <- "ATAC"
time10_cocl2_activity <- Signac::GeneActivity(time10_cocl2, 
                                       extend.downstream = 5000,
                                       extend.upstream = 5000,
                                       verbose = T)
set.seed(10)
Seurat::DefaultAssay(time10_dabtram) <- "ATAC"
time10_dabtram_activity <- Signac::GeneActivity(time10_dabtram, 
                                       extend.downstream = 5000,
                                       extend.upstream = 5000,
                                       verbose = T)

save(time0, time10_cis, time10_cocl2, time10_dabtram,
     time0_activity, time10_cis_activity,
     time10_cocl2_activity, time10_dabtram_activity,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_GEXATAC_and_ATACactivity.RData")

#################

time0b <- time0; Seurat::DefaultAssay(time0b) <- "RNA"; time0b[["ATAC"]] <- NULL
time10b_cis <- time10_cis; Seurat::DefaultAssay(time10b_cis) <- "RNA";  time10b_cis[["ATAC"]] <- NULL
time10b_cocl2 <- time10_cocl2; Seurat::DefaultAssay(time10b_cocl2) <- "RNA"; time10b_cocl2[["ATAC"]] <- NULL
time10b_dabtram <- time10_dabtram; Seurat::DefaultAssay(time10b_dabtram) <- "RNA"; time10b_dabtram[["ATAC"]] <- NULL

all_data <- merge(time0b, y = c(time10b_cis, time10b_cocl2, time10b_dabtram), 
                  add.cell.ids = c("time0", "time10_cis", "time10_cocl2", "time10_dabtram"), 
                  project = "All_Data", merge.data = T)

time0_atac <- Seurat::CreateSeuratObject(counts = time0_activity)
time10_cis_atac <- Seurat::CreateSeuratObject(counts = time10_cis_activity)
time10_cocl2_atac <- Seurat::CreateSeuratObject(counts = time10_cocl2_activity)
time10_dabtram_atac <- Seurat::CreateSeuratObject(counts = time10_dabtram_activity)

all_data_atac <- merge(time0_atac, y = c(time10_cis_atac, time10_cocl2_atac, time10_dabtram_atac), 
                       add.cell.ids = c("time0", "time10_cis", "time10_cocl2", "time10_dabtram"), 
                       project = "All_Data", merge.data = T)

all_data[["ATAC"]] <- all_data_atac[["RNA"]]
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

save(all_data, 
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_withATACactivity.RData")

################################
################################
################################

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


Seurat::DefaultAssay(all_data) <- "ATAC"
set.seed(10)
all_data <- Seurat::FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 500)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- intersect(jackpot_genes, rownames(all_data))
all_data[["ATAC"]]@var.features <- unique(c(jackpot_genes, all_data[["RNA"]]@var.features))
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE) 
all_data <- Seurat::RunUMAP(all_data, dims = 1:25)

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E")),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_atac-activity_umap_jackpot1.png"),
                plot1, device = "png", width = 12, height = 8, units = "in")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(setdiff(jackpot_genes, c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E"))),
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_atac-activity_umap_jackpot2.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")


plot1 <-Seurat::DimPlot(all_data, reduction = "umap",
                        group.by = "original_dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC activity score (LogNormalized),\n", length(all_data[["RNA"]]@var.features), " genes, using 30 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_alldata_atac-activity_500genes_umap_original_dataset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

