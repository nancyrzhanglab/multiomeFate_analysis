rm(list=ls())
library(Seurat)
library(Signac)
load("../../../../out/kevin/Writeup6h/Writeup6h_tcca_RNA-ATAC_DABTRAM.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

plot1 <- Seurat::DimPlot(all_data, reduction = "common_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_tcca-RNA-ATAC_DABTRAM_umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(all_data, reduction = "distinct1_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nRNA (Saver) distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_tcca-RNA-ATAC_DABTRAM_umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(all_data, reduction = "distinct2_tcca",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Tilted-CCA: RNA (Saver)+ATAC\nATAC distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_tcca-RNA-ATAC_DABTRAM_umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

######################

# do a UMAP of just the day10 cells
keep_vec <- rep(NA, ncol(all_data))
keep_vec[which(all_data$dataset == "day10_DABTRAM")] <- TRUE
all_data$keep <- keep_vec
all_data2 <- subset(all_data, keep == TRUE)

# multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
# dimred_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_mat")
# multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
# dimred_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_dimred")
# dimred <- cbind(dimred_1, dimred_2)
dimred <- multiSVD_obj[["tcca_obj"]]$common_score
dimred <- dimred[colnames(all_data2),]
set.seed(10)
seurat_umap <- Seurat::RunUMAP(dimred, 
                               assay = "RNA")
all_data2[["common_tcca"]] <- seurat_umap

treatment <- "DABTRAM"
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]

tier1_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier1_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier1_names <- colnames(all_data)[tier1_idx]
tier2_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier2_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier2_names <- colnames(all_data)[tier2_idx]
tier3_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% tier3_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == paste0("day10_", treatment))
)
tier3_names <- colnames(all_data)[tier3_idx]
tier_vec <- rep(NA, ncol(all_data))
names(tier_vec) <- colnames(all_data)
tier_vec[tier3_names] <- paste0("3high_winner_", treatment)
tier_vec[tier2_names] <- paste0("2mid_winner_", treatment)
tier_vec[tier1_names] <- paste0("1loser_", treatment)
all_data2$tier_vec <- tier_vec[colnames(all_data2)]

col_palette <- c("gray", "blue", "red")
names(col_palette) <- c("1loser_DABTRAM", "2mid_winner_DABTRAM", "3high_winner_DABTRAM")
names(col_palette) <- sort(unique(tier_vec))
p1 <- Seurat::DimPlot(all_data2, reduction = "common_tcca",
                      cols = col_palette,
                      group.by = "tier_vec")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_tcca_RNA-ATAC_DABTRAM_only-day-10_common.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

