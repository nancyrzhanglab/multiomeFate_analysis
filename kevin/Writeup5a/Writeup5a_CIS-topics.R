rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

keep_vec <- rep(0, ncol(all_data))
keep_vec[all_data$dataset %in% c("day0", "day10_CIS", "week5_CIS")] <- 1
all_data$keep <- keep_vec
all_data_subset <- subset(all_data, keep == 1)

feature_vec <- colnames(all_data_subset[["fasttopic_CIS"]]@cell.embeddings)
Seurat::Idents(all_data_subset) <- "dataset"
plot1 <- Seurat::VlnPlot(all_data_subset, 
                         features = feature_vec,
                         ncol = 6,
                         pt.size = 0)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup5a/Writeup5a_fasttopicsCIS-violin.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

######################


gene_list <- list(cellcycle =  c(cc.genes$s.genes, cc.genes$g2m.genes),
                  jackpot = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                              "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                              "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                              "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                              "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")))


# Custom heatmap
col_vec <- viridis::viridis(25)
break_vec <- seq(0, 1, length.out = 26)
for(i in 1:2){
  mat <- all_data_subset[["fasttopic_CIS"]]@feature.loadings
  idx <- which(rownames(mat) %in% gene_list[[i]])
  mat <- mat[idx,]
  rowname_vec <- rownames(mat)
  l1_vec <- apply(mat, 1, sum)
  mat <- diag(1/l1_vec) %*% mat
  rownames(mat) <- rowname_vec
  
  hclust_res <- stats::hclust(stats::dist(mat))
  mat <- mat[hclust_res$order,]
  
  ratio <- min(2000*(nrow(mat)-1)/(ncol(mat)-1), 6500)/(2000)
  png(paste0("../../../../out/figures/Writeup5a/Writeup5a_saver-CIS_fasttopics-scores_", names(gene_list)[i], ".png"),
      height = 2000*ratio, width = 2000, res = 500, units = "px")
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  image(tiltedCCA:::.rotate(mat), 
        asp = ratio, ylim = c(0,1),
        main = "", xlab = "", ylab = "",
        xaxt = "n", yaxt = "n", bty = "n",
        breaks = break_vec, col = col_vec)
  
  # draw lines
  x_indent_val <- 1/(2*(ncol(mat)-1))
  x_vec <- seq(-x_indent_val, 1+x_indent_val, by = 2*x_indent_val)
  for(x in x_vec[seq(31,6,by=-5)[-1]]){
    graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 2, col = "white")
    graphics::lines(y = c(0,1), x = rep(x, 2), lwd = 1.5, lty = 2)
  }
  
  # label genes
  y_indent_val <- 1/(2*(nrow(mat)-1))
  y_vec <- seq(-y_indent_val, 1+y_indent_val, by = 2*y_indent_val)
  for(y in y_vec[seq(6,length(y_vec), by = 5)]){
    graphics::lines(x = c(0,1), y = rep(y, 2), lwd = 2, col = "white")
  }
  x_vec_label <- rep(x_vec[c(3,5)], times = ceiling(nrow(mat)/2))[1:nrow(mat)]
  y_vec_label <- seq(1, 0, by = -2*y_indent_val)
  graphics::text(x = x_vec_label, y = y_vec_label, labels = rownames(mat),
                 col = "white", cex = 0.5)
  
  graphics.off()
}
