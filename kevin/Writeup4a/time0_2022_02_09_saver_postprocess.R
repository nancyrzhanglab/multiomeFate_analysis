rm(list=ls())
load("../../../../out/kevin/Writeup4a/time0_2022_02_09_saver.RData")
library(Seurat)
library(SAVER)

jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <-  intersect(rownames(saver_res$estimate), jackpot_genes)
round(sapply(jackpot_genes, function(gene_name){
  j <- which(rownames(saver_res$estimate) == gene_name)
  quantile(saver_res$estimate[j,])
}), 2)

cor_mat <- SAVER::cor.genes(saver_res)
jackpot_idx <- which(rownames(saver_res$estimate) %in% jackpot_genes)
round(cor_mat[jackpot_idx, jackpot_idx],2)
cor_mat_subset <- cor_mat[jackpot_idx, jackpot_idx]

diag(cor_mat_subset) <- 0
png("../../../../out/figures/Writeup4a/2022-02-09_saver_correlation.png",
    height = 1500, width = 1500, units = "px", res = 300)
gplots::heatmap.2(cor_mat_subset, scale = "none", col = gplots::bluered(100), 
          trace = "none", density.info = "none")
graphics.off()

t0_obj[["saver"]] <- Seurat::CreateAssayObject(counts = saver_res$estimate)
Seurat::DefaultAssay(t0_obj) <- "saver"
plot1 <- Seurat::FeatureScatter(t0_obj, 
                                feature1 = "FN1", feature2 = "EGFR",
                                group.by = NULL)
plot1 <- plot1 + ggplot2::theme_gray() + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_saver_fn1_egfr.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

Seurat::DefaultAssay(t0_obj) <- "saver"
plot1 <- Seurat::FeatureScatter(t0_obj, 
                                feature1 = "AXL", feature2 = "NDRG1",
                                group.by = NULL)
plot1 <- plot1 + ggplot2::theme_gray() + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_saver_axl_ndgr1.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

####################

diag(cor_mat) <- 0
idx <- which.max(cor_mat["FN1",])
feature2 <- rownames(cor_mat)[idx]
Seurat::DefaultAssay(t0_obj) <- "saver"
plot1 <- Seurat::FeatureScatter(t0_obj, 
                                feature1 = "FN1", feature2 = feature2,
                                group.by = NULL)
plot1 <- plot1 + ggplot2::theme_gray() + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_saver_fn1_mostcorrelated.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")


