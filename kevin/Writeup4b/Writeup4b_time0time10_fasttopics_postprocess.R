rm(list=ls())
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_fasttopics_onlyGEX.RData")

library(Seurat)
library(fastTopics)

fasttopic_dimred <- Seurat::CreateDimReducObject(topic_res$L, 
                                                 assay = "RNA",
                                                 key = "fastTopic_")
all_data[["fasttopic"]] <- fasttopic_dimred

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = paste0("fastTopic_", 1:ncol(topic_res$L)),
                             ncol = 6)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_time0time10_umap_fasttopics.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

#################

coef_mat <- topic_res$F
coef_mat <- scale(coef_mat, center = T, scale = T)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- sort(unique(intersect(jackpot_genes, rownames(coef_mat))))
round(apply(coef_mat[jackpot_genes,], 2, stats::quantile),2)
round(apply(coef_mat, 2, stats::quantile),2)
round(coef_mat[jackpot_genes,11],2)
round(coef_mat[jackpot_genes,])

round(quantile(coef_mat[,11], probs = c(0,0.5,0.95,0.99,0.995,1)),2)
rownames(coef_mat)[which(coef_mat[,11] > 9)]

##################3

coef_mat <- topic_res$F
coef_mat <- scale(coef_mat, center = T, scale = T)
max_val <- apply(abs(coef_mat), 1, max)
jackpot_idx <- which(rownames(coef_mat) %in% jackpot_genes)
val <- c(max_val[-jackpot_idx], max_val[jackpot_idx])
perturbation <- min(diff(sort(unique(max_val))))
idx <- rank(max_val+stats::runif(length(max_val), min = 0, max = 2*perturbation))
idx <- c(idx[-jackpot_idx], idx[jackpot_idx])
labeling_vec <- rep(0, nrow(coef_mat))
labeling_vec[jackpot_idx] <- 1
labeling_vec <- c(labeling_vec[-jackpot_idx], labeling_vec[jackpot_idx])
labeling_vec <- as.factor(labeling_vec)
text_vec <- rownames(coef_mat)
text_vec <- c(text_vec[-jackpot_idx], text_vec[jackpot_idx])

set.seed(10)
df <- data.frame(Max_fasttopic_score = val, Gene_rank = idx, 
                 labeling = labeling_vec,
                 text = text_vec)

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = Gene_rank, y = Max_fasttopic_score)) 
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == 1), 
                                          ggplot2::aes(label = text, color = labeling),
                                          box.padding = ggplot2::unit(0.5, 'lines'),
                                          point.padding = ggplot2::unit(1.6, 'lines'),
                                          size = 3,
                                          nudge_y = 2)
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_time0time10_fasttopics_magnitudes.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")


