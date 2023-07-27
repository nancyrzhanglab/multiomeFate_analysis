rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")
all_data2$keep <- !is.na(all_data2$assigned_lineage)
all_data2 <- subset(all_data2, keep == TRUE)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_DABTRAM",
                            dims = 1:30)

# let's first plot the data as-is
col_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(col_vec) <- c("week5_DABTRAM", "day10_DABTRAM", "day0")
plot1 <- Seurat::DimPlot(all_data, reduction = "umap", cols = col_vec,
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("DABTRAM (Only cells with assigned lineages)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

###################
cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")


lineage_name <- "Lin24619"
cell_idx <- intersect(which(all_data$dataset == "week5_DABTRAM"),
                      which(all_data$assigned_lineage == lineage_name))
cell_names <- colnames(all_data)[cell_idx]

cell_idx2 <- intersect(which(all_data$dataset == "day10_DABTRAM"),
                       which(all_data$assigned_lineage == lineage_name))
cell_names2 <- colnames(all_data)[cell_idx2]

cell_idx3 <- intersect(which(all_data$dataset == "day0"),
                       which(all_data$assigned_lineage == lineage_name))
cell_names3 <- colnames(all_data)[cell_idx3]

p1 <- Seurat::DimPlot(all_data, reduction = "umap",
                      cells.highlight = list(week5 = cell_names, day10 = cell_names2, day0 = cell_names3),
                      cols.highlight = cols.highlight_vec)
p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                                   length(cell_idx2), " day10 cells and ", length(cell_idx3), " day0 cells)")) +
  ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6k/Writeup6k_DABTRAM_umap_lineage-", lineage_name, ".png"),
                p1, device = "png", width = 6, height = 5, units = "in")

