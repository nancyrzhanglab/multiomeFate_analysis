rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-chromAct_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::Idents(all_data) = "dataset"
dabtram_only <- subset(all_data, idents=c("day0", "day10_DABTRAM", "week5_DABTRAM"))
dabtram_only <- Seurat::RunUMAP(dabtram_only, reduction.key="ftUMAP", reduction = "fasttopic_DABTRAM", dims=c(1:30))

# How many cells need to be seen to count as "win"?
min_cell_count <- 50
earlierday_cond_str <- "day10_DABTRAM"
later_cond <- "week5_DABTRAM"

lintab <- as.data.frame.matrix(table(all_data$assigned_lineage, all_data$dataset)) # lineage by condition table

sel.lineages <- rownames(lintab)[which(lintab[,later_cond] >= min_cell_count)]  # lineages that survived in given condition
cells.1 <- which(dabtram_only$assigned_lineage %in% sel.lineages 
                 & dabtram_only$assigned_posterior > 0.5
                 & dabtram_only$dataset == earlierday_cond_str)

############

latercol <- "chartreuse4"

plot1 <- Seurat::DimPlot(dabtram_only, reduction="umap", label=TRUE, label.size=7) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::geom_point(
  data = data.frame(x = dabtram_only@reductions[["umap"]]@cell.embeddings[cells.1,1],
                    y = dabtram_only@reductions[["umap"]]@cell.embeddings[cells.1,2]),
  ggplot2::aes(x=x,y=y), col=latercol, size=3, shape=17)
plot1 <- plot1 + ggplot2::ggtitle(paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond))

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-chromAct_DABTRAM_frontrunner-", min_cell_count, ".png"),
                plot1, device = "png", width = 7, height = 6, units = "in")

############

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

quantile(alignment_vec)
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

###########

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-chromAct_DABTRAM_frontrunner-", min_cell_count, "_synchrony.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = dabtram_only@reductions[["umap"]]@cell.embeddings[,1],
     y = dabtram_only@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Synchrony"),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")

points(x = dabtram_only@reductions[["umap"]]@cell.embeddings[cells.1,1],
       y = dabtram_only@reductions[["umap"]]@cell.embeddings[cells.1,2],
       col = "white", pch = 16, cex = 3.5)
points(x = dabtram_only@reductions[["umap"]]@cell.embeddings[cells.1,1],
       y = dabtram_only@reductions[["umap"]]@cell.embeddings[cells.1,2],
       col = color_vec[cells.1], pch = 16, cex = 2)
graphics.off()

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-chromAct_DABTRAM_synchrony.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = dabtram_only@reductions[["umap"]]@cell.embeddings[,1],
     y = dabtram_only@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Synchrony"),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
graphics.off()

#############################

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+chromAct\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-chromAct_DABTRAM_umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(all_data, reduction = "distinct1_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+chromAct\nRNA (Saver) distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-chromAct_DABTRAM_umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(all_data, reduction = "distinct2_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+chromAct\nchromAct distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-chromAct_DABTRAM_umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")
