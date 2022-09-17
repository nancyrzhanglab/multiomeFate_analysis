rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_tcca_RNA-ATAC.RData")

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup5a/Writeup5a_tcca_RNA-ATAC_umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "dataset")
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup5a/Writeup5a_tcca_RNA-ATAC_umap_common-cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

###########

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
all_data$alignment <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(all_data$alignment), max(all_data$alignment), length.out = num_color)
color_vec <- sapply(all_data$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})
names(alignment_vec) <- colnames(all_data)

png(paste0("../../../../out/figures/Writeup5a/Writeup5a_tcca_RNA-ATAC_steadystate-full_cleaned.png"),
    height = 2000, width = 2000, units = "px", res = 500)
par(mar = c(0.5,0.5,0.5,0.5))
plot(x = all_data[["common_tcca"]]@cell.embeddings[,1],
     y = all_data[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = "",
     xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

Seurat::Idents(all_data) <- "dataset"
plot1 <- Seurat::VlnPlot(all_data, 
                         features = "alignment")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup5a/Writeup5a_steadystate-violin.png"),
                plot1, device = "png", width = 10, height = 8, units = "in")


