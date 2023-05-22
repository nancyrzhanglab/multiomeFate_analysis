rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "COCL2"
relevant_celltypes <- c("day0", paste0(c("day10_", "week5_"), treatment))
all_data2 <- subset(all_data, dataset %in% relevant_celltypes)
Seurat::DefaultAssay(all_data2) <- "Saver"

all_data2[["geneActivity"]] <- NULL
all_data2[["fasttopic_CIS"]] <- NULL
all_data2[["fasttopic_DABTRAM"]] <- NULL
all_data2[["common_tcca"]] <- NULL
all_data2[["distinct1_tcca"]] <- NULL
all_data2[["distinct2_tcca"]] <- NULL

all_data2[["umap"]] <- NULL
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, 
                             reduction = "fasttopic_COCL2",
                             dims = 1:30)

# let's first plot the data as-is
col_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(col_vec) <- c("week5_COCL2", "day10_COCL2", "day0")
plot1 <- Seurat::DimPlot(all_data2, reduction = "umap", cols = col_vec,
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_COCL2_umap-as-is.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot week5 cells by their lineage
tab_mat <- table(all_data2$assigned_lineage, all_data2$dataset)
length(which(tab_mat[,"week5_COCL2"] >= 20))
cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")

tab_mat <- tab_mat[order(tab_mat[,"week5_COCL2"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_COCL2"] >= 20)]
pdf("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_COCL2_large-lineages_week5.pdf", 
    onefile = T, width = 5, height = 5)
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  
  lineage_name <- lineage_vec[i]
  cell_idx <- intersect(which(all_data2$dataset == paste0("week5_", treatment)),
                        which(all_data2$assigned_lineage == lineage_name))
  cell_idx <- intersect(cell_idx, 
                        which(all_data2$assigned_posterior >= 0.5))
  cell_names <- colnames(all_data2)[cell_idx]
  
  cell_idx2 <- intersect(which(all_data2$dataset == paste0("day10_", treatment)),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx2 <- intersect(cell_idx2, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names2 <- colnames(all_data2)[cell_idx2]
  
  cell_idx3 <- intersect(which(all_data2$dataset == "day0"),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx3 <- intersect(cell_idx3, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names3 <- colnames(all_data2)[cell_idx3]
  
  p1 <- Seurat::DimPlot(all_data2, reduction = "umap",
                        cells.highlight = list(week5 = cell_names, day10 = cell_names2, day0 = cell_names3),
                        cols.highlight = cols.highlight_vec)
  p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                                     length(cell_idx2), " day10 cells and ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()

# do the same for day10
tab_mat <- tab_mat[order(tab_mat[,"day10_COCL2"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"day10_COCL2"] >= 20)]
pdf("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_COCL2_large-lineages_day10.pdf", 
    onefile = T, width = 5, height = 5)
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  
  lineage_name <- lineage_vec[i]
  cell_idx <- intersect(which(all_data2$dataset == paste0("week5_", treatment)),
                        which(all_data2$assigned_lineage == lineage_name))
  cell_idx <- intersect(cell_idx, 
                        which(all_data2$assigned_posterior >= 0.5))
  cell_names <- colnames(all_data2)[cell_idx]
  
  cell_idx2 <- intersect(which(all_data2$dataset == paste0("day10_", treatment)),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx2 <- intersect(cell_idx2, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names2 <- colnames(all_data2)[cell_idx2]
  
  cell_idx3 <- intersect(which(all_data2$dataset == "day0"),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx3 <- intersect(cell_idx3, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names3 <- colnames(all_data2)[cell_idx3]
  
  if(length(cell_names) == 0){
    p1 <- Seurat::DimPlot(all_data2, reduction = "umap",
                          cells.highlight = list(day10 = cell_names2, day0 = cell_names3),
                          cols.highlight = cols.highlight_vec[c("day10", "day0")])
  } else {
    p1 <- Seurat::DimPlot(all_data2, reduction = "umap",
                          cells.highlight = list(week5 = cell_names, day10 = cell_names2, day0 = cell_names3),
                          cols.highlight = cols.highlight_vec)
  }
  p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: ", lineage_name, " (", length(cell_idx2), " day10 cells,\nto ",
                                     length(cell_idx), " week5 cells from ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()

####################################

low_posterior_cells <- colnames(all_data2)[which(all_data2$assigned_posterior <= 0.5)]

# let's first plot the data as-is

low_posterior_cells_subset <- intersect(low_posterior_cells,
                                        colnames(all_data2)[which(all_data2$dataset == "day0")])
p1 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                         cells.highlight = low_posterior_cells_subset)
p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: Low posterior (Day 0)"))
p1 <- p1 + Seurat::NoLegend()

low_posterior_cells_subset <- intersect(low_posterior_cells,
                                        colnames(all_data2)[which(all_data2$dataset == "day10_COCL2")])
p2 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = low_posterior_cells_subset)
p2 <- p2 + ggplot2::ggtitle(paste0("COCL2: Low posterior (Day 10)"))
p2 <- p2 + Seurat::NoLegend()

low_posterior_cells_subset <- intersect(low_posterior_cells,
                                        colnames(all_data2)[which(all_data2$dataset == "week5_COCL2")])
p3 <- Seurat::DimPlot(all_data2, reduction = "umap", 
                      cells.highlight = low_posterior_cells_subset)
p3 <- p3 + ggplot2::ggtitle(paste0("COCL2: Low posterior (Week 5)"))
p3 <- p3 + Seurat::NoLegend()

p <- cowplot::plot_grid(p1, p2, p3, ncol = 3)

ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_COCL2_umap_low-posterior.png"),
                p, device = "png", width = 15, height = 5, units = "in")

###################

all_data3 <- subset(all_data2, dataset == "day10_COCL2")

treatment <- "COCL2"
# find the winning and losing cells
tab_mat <- table(all_data2$assigned_lineage, all_data2$dataset)
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 3),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 49))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]

tier1_idx <- intersect(
  intersect(which(all_data3$assigned_lineage %in% tier1_lineages),
            which(all_data3$assigned_posterior >= 0.5)),
  which(all_data3$dataset == paste0("day10_", treatment))
)
tier2_idx <- intersect(
  intersect(which(all_data3$assigned_lineage %in% tier2_lineages),
            which(all_data3$assigned_posterior >= 0.5)),
  which(all_data3$dataset == paste0("day10_", treatment))
)
tier3_idx <- intersect(
  intersect(which(all_data3$assigned_lineage %in% tier3_lineages),
            which(all_data3$assigned_posterior >= 0.5)),
  which(all_data3$dataset == paste0("day10_", treatment))
)
length(tier1_idx); length(tier2_idx); length(tier3_idx)

tier_vec <- rep(NA, ncol(all_data3))
tier_vec[tier1_idx] <- paste0("1loser_", treatment)
tier_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
tier_vec[tier3_idx] <- paste0("3high_winner_", treatment)
table(tier_vec); table(is.na(tier_vec))
all_data3$tier <- tier_vec

set.seed(10)
all_data3 <- Seurat::RunPCA(all_data3, assay = "Saver", verbose = F)
set.seed(10)
all_data3 <- Seurat::RunUMAP(all_data3, 
                             reduction = "pca",
                             dims = 1:30)

col_palette <- c("gray", "blue", "red")
names(col_palette) <- sort(unique(all_data3$tier))
p1 <- Seurat::DimPlot(all_data3, reduction = "umap",
                      cols = col_palette,
                      group.by = "tier")
p1 <- p1 + ggplot2::ggtitle("COCL2 Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_COCL2_umap_day10.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

###################################################

col_palette <- c("gray", "blue", "red")
names(col_palette) <- sort(unique(all_data3$tier))
cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")

# do the same for day10
tab_mat <- tab_mat[order(tab_mat[,"week5_COCL2"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_COCL2"] >= 20)]
pdf("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_COCL2_week5_withday10umap.pdf", 
    onefile = T, width = 10, height = 5)
for(i in 1:length(lineage_vec)){
  print(paste0(i, " out of ", length(lineage_vec)))
  
  lineage_name <- lineage_vec[i]
  cell_idx <- intersect(which(all_data2$dataset == paste0("week5_", treatment)),
                        which(all_data2$assigned_lineage == lineage_name))
  cell_idx <- intersect(cell_idx, 
                        which(all_data2$assigned_posterior >= 0.5))
  cell_names <- colnames(all_data2)[cell_idx]
  
  cell_idx2 <- intersect(which(all_data2$dataset == paste0("day10_", treatment)),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx2 <- intersect(cell_idx2, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names2 <- colnames(all_data2)[cell_idx2]
  
  cell_idx3 <- intersect(which(all_data2$dataset == "day0"),
                         which(all_data2$assigned_lineage == lineage_name))
  cell_idx3 <- intersect(cell_idx3, 
                         which(all_data2$assigned_posterior >= 0.5))
  cell_names3 <- colnames(all_data2)[cell_idx3]
  
  if(length(cell_names) == 0){
    p1 <- Seurat::DimPlot(all_data2, reduction = "umap",
                          cells.highlight = list(day10 = cell_names2, day0 = cell_names3),
                          cols.highlight = cols.highlight_vec[c("day10", "day0")])
  } else {
    p1 <- Seurat::DimPlot(all_data2, reduction = "umap",
                          cells.highlight = list(week5 = cell_names, day10 = cell_names2, day0 = cell_names3),
                          cols.highlight = cols.highlight_vec)
  }
  p1 <- p1 + ggplot2::ggtitle(paste0("COCL2: ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                                     length(cell_idx2), " day10 cells and ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
 
  p2 <- Seurat::DimPlot(all_data3, reduction = "umap",
                        cells.highlight = cell_names2)
  p2 <- p2 + ggplot2::ggtitle("COCL2 UMAP among Day10s") + Seurat::NoLegend()
  
  p <- cowplot::plot_grid(p1, p2, ncol = 2)
  print(p)
}

dev.off()


