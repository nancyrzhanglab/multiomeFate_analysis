rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-chromAct_DABTRAM.RData")

# find the winning and losing cells
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
surviving_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_DABTRAM")] >= 50)]
dying_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_DABTRAM")] <= 1)]
winning_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% surviving_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day10_DABTRAM")
)
dying_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% dying_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day10_DABTRAM")
)
length(winning_idx); length(dying_idx)

# set idents
keep_vec <- rep(NA, ncol(all_data))
keep_vec[winning_idx] <- "winning_week5"
keep_vec[dying_idx] <- "dying_week5"
all_data$keep <- keep_vec
Seurat::Idents(all_data) <- "keep"

Seurat::DefaultAssay(all_data) <- "geneActivity"
gene_vec <- sort(intersect(rownames(multiSVD_obj$svd_1$v), rownames(multiSVD_obj$svd_2$v)))
de_res <- Seurat::FindMarkers(all_data,
                              features = gene_vec,
                              ident.1 = "winning_week5",
                              ident.2 = "dying_week5",
                              test.use = "wilcox",
                              slot = "scale.data",
                              min.pct = 0,
                              logfc.threshold = 0,
                              only.pos = F,
                              verbose = T)
logpval_vec <- -log10(de_res[,"p_val"])
names(logpval_vec) <- rownames(de_res)

##################

rsquare_vec <- tiltedCCA:::postprocess_modality_alignment(input_obj = multiSVD_obj,
                                                          bool_use_denoised = T,
                                                          seurat_obj = all_data,
                                                          input_assay = 1,
                                                          seurat_assay = "geneActivity",
                                                          seurat_slot = "data")

zz <- sort(intersect(names(rsquare_vec), names(logpval_vec)))
logpval_vec <- logpval_vec[zz]
rsquare_vec <- rsquare_vec[zz]

##################

source("../Writeup6b/gene_list.R")
important_genes <- sort(keygenes$DABTRAM)
important_genes <- intersect(important_genes, gene_vec)
length(important_genes)

png("../../../../out/figures/Writeup6d/Writeup6d_DABTRAM_RNA-ChromAct_differential-alignment.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "DABTRAM Week5 (RNA+ChromAct)\nGene differentiability vs. alignment",
                           bool_mark_ymedian = T,
                           col_gene_highlight_border = rgb(255, 205, 87, 255*0.5, maxColorValue = 255),
                           col_points = rgb(0.6, 0.6, 0.6, 0.1),
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           lty_polygon = 2,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon_bold = 5,
                           mark_median_xthres = 0)
graphics.off()

png("../../../../out/figures/Writeup6d/Writeup6d_DABTRAM_RNA-ChromAct_differential-alignment_important-genes.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "DABTRAM Week5 (RNA+ChromAct)\nGene differentiability vs. alignment",
                           bool_mark_ymedian = F,
                           bool_polygon_mean = T,
                           col_points = rgb(0.5, 0.5, 0.5, 0.1),
                           col_gene_highlight = 2,
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           density = 10,
                           gene_names = important_genes,
                           lty_polygon = 1,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon = 2,
                           lwd_polygon_bold = 4)
graphics.off()



