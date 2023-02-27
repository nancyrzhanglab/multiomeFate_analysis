rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_COCL2.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::Idents(all_data) = "dataset"
all_data <- Seurat::RunUMAP(all_data, reduction.key="ftUMAP", reduction = "fasttopic_COCL2", dims=c(1:30))

# How many cells need to be seen to count as "win"?
min_cell_count <- 50
earlierday_cond_str <- "day10_COCL2"
later_cond <- "week5_COCL2"

lintab <- as.data.frame.matrix(table(all_data$assigned_lineage, all_data$dataset)) # lineage by condition table

sel.lineages <- rownames(lintab)[which(lintab[,later_cond] >= min_cell_count)]  # lineages that survived in given condition
cells.1 <- which(all_data$assigned_lineage %in% sel.lineages 
                 & all_data$assigned_posterior > 0.5
                 & all_data$dataset == earlierday_cond_str)

############

latercol <- "chartreuse4"

plot1 <- Seurat::DimPlot(all_data, reduction="umap", label=TRUE, label.size=7) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::geom_point(
  data = data.frame(x = all_data@reductions[["umap"]]@cell.embeddings[cells.1,1],
                    y = all_data@reductions[["umap"]]@cell.embeddings[cells.1,2]),
  ggplot2::aes(x=x,y=y), col=latercol, size=3, shape=17)
plot1 <- plot1 + ggplot2::ggtitle(paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond))

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_frontrunner-", min_cell_count, ".png"),
                plot1, device = "png", width = 9, height = 6, units = "in")

##############

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

##############

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_frontrunner-", min_cell_count, "_synchrony.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = all_data@reductions[["umap"]]@cell.embeddings[,1],
     y = all_data@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Synchrony"),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "",
     cex.main = 0.5)

points(x = all_data@reductions[["umap"]]@cell.embeddings[cells.1,1],
       y = all_data@reductions[["umap"]]@cell.embeddings[cells.1,2],
       col = "white", pch = 16, cex = 3.5)
points(x = all_data@reductions[["umap"]]@cell.embeddings[cells.1,1],
       y = all_data@reductions[["umap"]]@cell.embeddings[cells.1,2],
       col = color_vec[cells.1], pch = 16, cex = 2)
graphics.off()

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_synchrony.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(x = all_data@reductions[["umap"]]@cell.embeddings[,1],
     y = all_data@reductions[["umap"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = paste("Cells in",earlierday_cond_str,"whose lineage has >",min_cell_count,"cells in",later_cond,": Synchrony"),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "",
     cex.main = 0.5)
graphics.off()

cells.2 <- which(!all_data$assigned_lineage %in% sel.lineages 
                 & all_data$assigned_posterior > 0.5
                 & all_data$dataset == earlierday_cond_str)

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_synchrony-histogram.png"),
    height = 1500, width = 3000, units = "px", res = 500)
par(mfrow = c(1,2), mar = c(4,4,4,0.5))

hist(all_data$alignment[cells.1], breaks = 50, col = "gray", 
     main = paste0("Cells in ",earlierday_cond_str," whose lineage\nhas > ",min_cell_count," cells in ",later_cond),
     cex.main = 0.6,
     xlim = range(all_data$alignment))
mean_val <- mean(all_data$alignment[cells.1])
median_val <- median(all_data$alignment[cells.1])
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2, lty = 1)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

hist(all_data$alignment[cells.2], breaks = 50, col = "gray", 
     main = paste0("Cells in ",earlierday_cond_str," whose lineage\nhas < ",min_cell_count," cells in ",later_cond),
     cex.main = 0.6,
     xlim = range(all_data$alignment))
mean_val <- mean(all_data$alignment[cells.2])
median_val <- median(all_data$alignment[cells.2])
lines(rep(mean_val, 2), c(-1e6,1e6), col = 2, lwd = 2, lty = 1)
lines(rep(median_val, 2), c(-1e6,1e6), col = 3, lwd = 2, lty = 2)

graphics.off()

set.seed(10)
t_test <- stats::t.test(x = all_data$alignment[cells.1],
                        y = all_data$alignment[cells.2],
                        alternative = "two.sided")
wilcox_test <- stats::wilcox.test(x = all_data$alignment[cells.1],
                                  y = all_data$alignment[cells.2],
                                  alternative = "two.sided")
# COCL2
t_test
wilcox_test

#################################################

# make a plot of alignment-score against number-of-cells-in-future

cell_idx <- intersect(intersect(which(all_data$dataset == "day10_COCL2"),
                                which(all_data$assigned_posterior > 0.5)),
                      which(!is.na(all_data$assigned_lineage)))
alignment_vec <- all_data$alignment[cell_idx]
posterior_vec <- (all_data$assigned_posterior > 0.5)[cell_idx]
lineage_vec <- all_data$assigned_lineage[cell_idx]

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
future_num_vec <- sapply(lineage_vec, function(lineage){
  tab_mat[lineage,"week5_COCL2"]
})

number_uniq <- sort(unique(future_num_vec))
average_synchrony <- sapply(number_uniq, function(val){
  median(alignment_vec[which(future_num_vec == val)])
})

png(paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_synchrony-vs-number.png"),
    height = 2000, width = 2000, units = "px", res = 500)
plot(alignment_vec, log1p(future_num_vec),
     xlab = "COCL2 day10 synchrony",
     ylab = "Log1p # week5 cells in lineage",
     pch = 16, col = rgb(0.5,0.5,0.5,0.2))

for(i in 1:length(number_uniq)){
  points(average_synchrony[i], log1p(number_uniq[i]),
         col = "white", cex = 2, pch = 16)
  points(average_synchrony[i], log1p(number_uniq[i]),
         col = 2, cex = 1.5, pch = 16)
}

lines(c(-1e2,1e2), rep(log1p(50),2), col = 3, lwd = 2, lty = 2)
graphics.off()

#############################

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(all_data, reduction = "distinct1_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nRNA (Saver) distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(all_data, reduction = "distinct2_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nATAC distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca-RNA-ATAC_COCL2_umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")


