rm(list=ls())
load("../../../../out/kevin/Writeup4a/2022-02-09_time0_preprocessed.RData")
library(Seurat)
library(SAVER)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <-  intersect(rownames(t0_obj[["RNA"]]), jackpot_genes)
gene_vec <- unique(c(Seurat::VariableFeatures(t0_obj, assay = "SCT"),
                     jackpot_genes))
mat <- t0_obj[["RNA"]]@counts[Seurat::VariableFeatures(t0_obj, assay = "SCT"),]
saver_res <- SAVER::saver(x = mat, ncores = 1)
save(saver_res, t0_obj, date_of_run, session_info, gene_vec,
     file = "../../../../out/kevin/Writeup4a/time0_2022_02_09_saver.RData")

