rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day10_lineage-imputation_stepwise-up_step2.RData")
load("../../../../out/kevin/Writeup6i/Writeup6i_COCL2-day10_extracted.RData")
all_data2$tier_vec <- all_data2$keep

###################

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_COCL2",
                            dims = 1:30)

cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")


lineage_name <- "Lin5217"
cell_idx <- intersect(which(all_data$dataset == "week5_COCL2"),
                      which(all_data$assigned_lineage == lineage_name))
cell_names <- colnames(all_data)[cell_idx]

cell_idx2 <- intersect(which(all_data$dataset == "day10_COCL2"),
                       which(all_data$assigned_lineage == lineage_name))
cell_names2 <- colnames(all_data)[cell_idx2]

cell_idx3 <- intersect(which(all_data$dataset == "day0"),
                       which(all_data$assigned_lineage == lineage_name))
cell_names3 <- colnames(all_data)[cell_idx3]

p1 <- Seurat::DimPlot(all_data, reduction = "umap",
                      cells.highlight = list(week5 = cell_names, day10 = cell_names2, day0 = cell_names3),
                      cols.highlight = cols.highlight_vec)
p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                                   length(cell_idx2), " day10 cells and ", length(cell_idx3), " day0 cells)")) +
  ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_COCL2_umap_lineage-", lineage_name, ".png"),
                p1, device = "png", width = 6, height = 5, units = "in")
