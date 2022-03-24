rm(list=ls())
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_processed2.RData")
library(Seurat)
library(SAVER)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data) <- "RNA"
all_data <- Seurat::FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 2000)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- intersect(jackpot_genes, rownames(all_data2))
all_data2[["RNA"]]@var.features <- unique(c(jackpot_genes, all_data2[["RNA"]]@var.features))

mat <- all_data[["RNA"]]@counts[Seurat::VariableFeatures(all_data, assay = "RNA"),]
