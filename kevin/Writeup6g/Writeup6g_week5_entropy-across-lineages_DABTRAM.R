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


treatment <- "DABTRAM"
relevant_celltypes <- c("day0", paste0(c("day10_", "week5_"), treatment))
all_data2 <- subset(all_data, dataset %in% relevant_celltypes)
Seurat::DefaultAssay(all_data2) <- "Saver"

all_data2[["geneActivity"]] <- NULL
all_data2[["fasttopic_CIS"]] <- NULL
all_data2[["fasttopic_COCL2"]] <- NULL
all_data2[["common_tcca"]] <- NULL
all_data2[["distinct1_tcca"]] <- NULL
all_data2[["distinct2_tcca"]] <- NULL
all_data2[["activityPCA"]] <- NULL

all_data2[["umap"]] <- NULL
set.seed(10)
all_data2 <- Seurat::RunUMAP(all_data2, 
                             reduction = "fasttopic_DABTRAM",
                             dims = 1:30)

# let's first plot the data as-is
col_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(col_vec) <- c("week5_DABTRAM", "day10_DABTRAM", "day0")
plot1 <- Seurat::DimPlot(all_data2, reduction = "umap", cols = col_vec,
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("DABTRAM"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_DABTRAM_umap-as-is.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot week5 cells by their lineage
tab_mat <- table(all_data2$assigned_lineage, all_data2$dataset)
length(which(tab_mat[,"week5_DABTRAM"] >= 20))
cols.highlight_vec <- c("darkblue", "#DE2D26", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10", "day0")

tab_mat <- tab_mat[order(tab_mat[,"week5_DABTRAM"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"week5_DABTRAM"] >= 20)]
pdf("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_DABTRAM_large-lineages_week5.pdf", 
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
  p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                                     length(cell_idx2), " day10 cells and ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()

# do the same for day10
tab_mat <- tab_mat[order(tab_mat[,"day10_DABTRAM"], decreasing = T),]
lineage_vec <- rownames(tab_mat)[which(tab_mat[,"day10_DABTRAM"] >= 20)]
pdf("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_DABTRAM_large-lineages_day10.pdf", 
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
    p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: ", lineage_name, " (", length(cell_idx2), " day10 cells,\nto ",
                                     length(cell_idx), " week5 cells from ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()


# now plot the forerunners
dimred <- all_data2[["umap"]]@cell.embeddings
forerunner_idx <- intersect(which(dimred[,1] >= 4),
                            which(dimred[,2] >= -2.5))
forerunner_idx <- intersect(forerunner_idx, which(all_data2$dataset == paste0("day10_", treatment)))
forerunner_idx <- intersect(forerunner_idx, which(all_data2$assigned_posterior >= 0.5))

num_lineage <- sapply(rownames(tab_mat), function(lineage_name){
  cell_idx <- which(all_data2$assigned_lineage == lineage_name)
  length(intersect(cell_idx, forerunner_idx))
})
num_lineage <- num_lineage[order(num_lineage, decreasing = T)]
lineage_vec <- names(num_lineage)[num_lineage >= 5]

pdf("../../../../out/figures/Writeup6g/Writeup6g_week5_entropy-across-lineage_DABTRAM_large-lineages_day10-forerunners.pdf", 
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
  p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: ", lineage_name, " (", length(cell_idx2), " day10 cells,\nto ",
                                     length(cell_idx), " week5 cells from ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()

###################

all_data3 <- subset(all_data2, dataset == "day10_DABTRAM")

treatment <- "DABTRAM"
# find the winning and losing cells
tab_mat <- table(all_data2$assigned_lineage, all_data2$dataset)
tier1_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 50)]
tier2_lineages <- rownames(tab_mat)[intersect(which(tab_mat[,paste0("week5_", treatment)] >= 5),
                                              which(tab_mat[,paste0("week5_", treatment)] <= 25))]
tier3_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] <= 2)]
length(tier1_lineages); length(tier2_lineages); length(tier3_lineages)

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

keep_vec <- rep(NA, ncol(all_data3))
keep_vec[tier1_idx] <- paste0("3high_winner_", treatment)
keep_vec[tier2_idx] <- paste0("2mid_winner_", treatment)
keep_vec[tier3_idx] <- paste0("1loser_", treatment)
table(keep_vec)
all_data3$keep <- keep_vec
all_data3 <- subset(all_data3, keep == c(paste0("3high_winner_", treatment), 
                                         paste0("2mid_winner_", treatment),
                                         paste0("1loser_", treatment)))


all_data3[["umap"]] <- NULL
set.seed(10)
all_data3 <- Seurat::RunUMAP(all_data3, 
                             reduction = "fasttopic_DABTRAM",
                             dims = 1:30)

col_palette <- c("gray", "blue", "red")
names(col_palette) <- sort(unique(all_data3$keep))
p1 <- Seurat::DimPlot(all_data3, reduction = "umap",
                      cols = col_palette,
                      group.by = "keep")
p1 <- p1 + ggplot2::ggtitle("DABTRAM Day 10") 
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6g/Writeup6g_DABTRAM_umap_day10.png"),
                p1, device = "png", width = 7, height = 5, units = "in")




