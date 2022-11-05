rm(list=ls())
library(Seurat)
library(Signac)
library(igraph)

load("../../../../out/kevin/Writeup5a/Writeup5a_tcca_RNA-geneActivity.RData")
source("../Writeup5a/color_palette.R")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tab_mat <- log10(tab_mat+1)

png(paste0("../../../../out/figures/Writeup6b/Writeup6b_lineage_pairs.png"),
    height = 3000, width = 3000, res = 500, units = "px")
# pairs(tab_mat, lower.panel = NULL, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
psych::pairs.panels(tab_mat, 
                    stars = T, 
                    ellipses = F,
                    smooth = F,
                    lm = T,
                    smoother = F,
                    show.points = T,
                    cex.cor = 1,
                    cex = 0.5,
                    pch = 16,
                    gap = 0,
                    col = rgb(0.5,0.5,0.5,0.2))
graphics.off()

##################

# violin plots
treatment_vec <- c("CIS", "COCL2", "DABTRAM")
Seurat::Idents(all_data) <- "dataset"

for(treatment_name in treatment_vec){
  print(treatment_name)
  all_data2 <- all_data
  keep_vec <- rep(0, ncol(all_data2))
  keep_vec[which(all_data2$dataset %in% c("day0", paste0(c("day10_", "week5_"), treatment_name)))] <- 1
  all_data2$keep <- keep_vec
  all_data2 <- subset(all_data2, keep == 1)
  
  plot1 <- Seurat::VlnPlot(all_data2,
                           features = paste0("fastTopic", treatment_name, "_", 1:ncol(all_data2[[paste0("fasttopic_", treatment_name)]]@cell.embeddings)),
                           ncol = 5,
                           pt.size = 0)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6b/Writeup6b_", treatment_name, "_fastTopics.png"),
                  plot1, device = "png", width = 15, height = 20, units = "in")
}

##################################

weight_cis <- all_data[["fasttopic_CIS"]]@feature.loadings
weight_cocl2 <- all_data[["fasttopic_COCL2"]]@feature.loadings
weight_cis <- weight_cis[sort(intersect(rownames(weight_cis), rownames(weight_cocl2))),]
weight_cocl2 <- weight_cocl2[sort(intersect(rownames(weight_cis), rownames(weight_cocl2))),]

# weight_cis <- scale(weight_cis)
# weight_cocl2 <- scale(weight_cocl2)
zz <- cor(weight_cis, weight_cocl2)


