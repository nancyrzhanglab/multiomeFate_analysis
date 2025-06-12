rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6o/Writeup6o_extract_all-data_atac-peakvi_exploration.RData")
lsi <- all_data[["lsi"]]@cell.embeddings
library_count <- all_data$nCount_ATAC

png("../../../../out/figures/kevin/Writeup6o/Writeup6o_lsi.png",
    height = 5000, width = 1500, units = "px", res = 300)
par(mfrow = c(12,1), mar = c(0.5,4,4,0.5))
for(i in 1:12){
  hist(lsi[,i], breaks = 50, col = "gray",
       main = paste0("Latent dimension ", i, ": Corr = ", round(stats::cor(library_count, lsi[,i]), 2)),
       xlab = "")
  quantile_vec <- stats::quantile(lsi[,i], probs = c(0.25,0.5,0.75))
  lines(x = rep(quantile_vec[2],2), y = c(0,1e5), col = 2, lwd = 1.5)
  lines(x = rep(quantile_vec[1],2), y = c(0,1e5), col = 2, lwd = 1, lty = 2)
  lines(x = rep(quantile_vec[3],2), y = c(0,1e5), col = 2, lwd = 1, lty = 2)
}
graphics.off()

#################################################

treatment_vec <- c("CIS", "COCL2", "DABTRAM")

for(treatment in treatment_vec){
  print(treatment)
  peakvi_mat <- read.csv(paste0("../../../../out/kevin/Writeup6o/Writeup6o_all-data-atac_", treatment, "_peakVI.csv"),
                         row.names = 1)
  svd_res <- svd(peakvi_mat)
  tmp <- sweep(svd_res$u, MARGIN = 2, STATS = svd_res$d, FUN = "*")
  rownames(tmp) <- rownames(peakvi_mat)
  peakvi_mat <- tmp
  peakvi_mat <- scale(peakvi_mat)
  library_count <- all_data$nCount_ATAC[rownames(peakvi_mat)]

  png(paste0("../../../../out/figures/kevin/Writeup6o/Writeup6o_peakvi-", treatment, ".png"),
      height = 5000, width = 1500, units = "px", res = 300)
  par(mfrow = c(12,1), mar = c(0.5,4,4,0.5))
  for(i in 1:12){
    hist(peakvi_mat[,i], breaks = 50, col = "gray",
         main = paste0("Latent dimension ", i, ": Corr = ", round(stats::cor(library_count, peakvi_mat[,i]), 2)),
         xlab = "")
    quantile_vec <- stats::quantile(peakvi_mat[,i], probs = c(0.25,0.5,0.75))
    lines(x = rep(quantile_vec[2],2), y = c(0,1e5), col = 2, lwd = 1.5)
    lines(x = rep(quantile_vec[1],2), y = c(0,1e5), col = 2, lwd = 1, lty = 2)
    lines(x = rep(quantile_vec[3],2), y = c(0,1e5), col = 2, lwd = 1, lty = 2)
  }
  graphics.off()

  all_data2 <- all_data
  keep_vec <- rep(FALSE, ncol(all_data2))
  names(keep_vec) <- colnames(all_data2)
  keep_vec[rownames(peakvi_mat)] <- TRUE
  all_data2$keep <- keep_vec
  all_data2 <- subset(all_data2, keep == TRUE)

  all_data2[["peakvi"]] <- Seurat::CreateDimReducObject(embeddings = peakvi_mat,
                                                        key = "peakvi",
                                                        assay = "Saver")
  set.seed(10)
  all_data2 <- Seurat::RunUMAP(all_data2,
                               dims = 1:ncol(peakvi_mat),
                               reduction = "peakvi",
                               reduction.name = "atacumap")

  plot1 <-Seurat::DimPlot(all_data2, reduction = "atacumap",
                          group.by = "dataset", label = TRUE,
                          repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (Peak-VI), ", treatment))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/kevin/Writeup6o/Writeup6o_peakvi_inSeurat_umap_", treatment, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")

}
