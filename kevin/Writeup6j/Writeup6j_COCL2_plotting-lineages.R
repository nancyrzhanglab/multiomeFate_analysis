rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment <- "COCL2"

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

keep_vec <- rep(NA, ncol(all_data))
idx1 <- which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))
idx2 <- which(all_data$assigned_posterior >= 0.5)
idx3 <- which(!is.na(all_data$assigned_lineage))
keep_vec[intersect(intersect(idx1, idx2), idx3)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

all_data[["umap"]] <- NULL
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = "fasttopic_COCL2",
                            dims = 1:30)

# let's first plot the data as-is
col_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(col_vec) <- c("week5_COCL2", "day10_COCL2", "day0")
plot1 <- Seurat::DimPlot(all_data, reduction = "umap", cols = col_vec,
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6j/Writeup6j_COCL2_umap_only-confident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#############

# plot week5 cells by their lineage
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
length(which(tab_mat[,"week5_COCL2"] >= 20))
cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")

tab_mat <- tab_mat[order(tab_mat[,"week5_COCL2"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_COCL2"] >= 10)]
pdf("../../../../out/figures/Writeup6j/Writeup6j_COCL2-lineages_week5.pdf", 
    onefile = T, width = 5, height = 5)
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  
  lineage_name <- lineage_vec[i]
  cell_idx <- intersect(which(all_data$dataset == paste0("week5_", treatment)),
                        which(all_data$assigned_lineage == lineage_name))
  cell_names <- colnames(all_data)[cell_idx]
  
  cell_idx2 <- intersect(which(all_data$dataset == paste0("day10_", treatment)),
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
  print(p1)
}

dev.off()