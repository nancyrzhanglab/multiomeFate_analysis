rm(list=ls())
load("../../../../out/kevin/Writeup4a/time0_2022_02_09_fasttopics.RData")
load("../../../../out/kevin/Writeup4a/2022-02-09_time0_preprocessed.RData")
library(Seurat)
library(fastTopics)

fasttopic_dimred <- Seurat::CreateDimReducObject(topic_res$L, 
                                                 assay = "SCT",
                                                 key = "fastTopic_")
t0_obj[["fasttopic"]] <- fasttopic_dimred

plot1 <- Seurat::FeaturePlot(t0_obj, 
                             reduction = "umap",
                             features = paste0("fastTopic_", 1:ncol(topic_res$L)),
                             ncol = 6)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_umap_fasttopics.png"),
                plot1, device = "png", width = 20, height = 10, units = "in")

#################

coef_mat <- topic_res$F
coef_mat <- scale(coef_mat, center = T, scale = T)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- unique(intersect(jackpot_genes, rownames(coef_mat)))
round(apply(coef_mat[jackpot_genes,], 2, stats::quantile),2)
round(apply(coef_mat, 2, stats::quantile),2)
round(coef_mat[jackpot_genes,2],2)

