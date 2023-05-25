rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")

cell_names <- colnames(all_data2)
tier_vec <- all_data2$tier_vec[cell_names]
tier_vec <- tier_vec[!is.na(tier_vec)]
cell_names <- names(tier_vec)
pseudotime <- all_data$pseudotime[cell_names]

# let's first make a violin plot
df <- data.frame(tier = tier_vec,
                 pseudotime = pseudotime)
col_vec <- c("lightgray", "#7190AF", "dodgerblue4")
names(col_vec) <- c("1loser_DABTRAM", "2mid_winner_DABTRAM", "3high_winner_DABTRAM")

p1 <- ggplot2::ggplot(df, ggplot2::aes(x=tier, y=pseudotime)) +
  ggplot2::geom_violin(trim=FALSE, ggplot2::aes(fill=tier)) +
  ggplot2::scale_fill_manual(values=col_vec) +
  ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::geom_boxplot(width=0.05)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime_violinplot.png"),
                p1, device = "png", width = 7, height = 5, units = "in")

# where are the day0 cells that do not expand, and how many unique lineages do they belong to?
day0_cells <- intersect(names(tier_vec)[tier_vec == "1loser_DABTRAM"],
                        names(pseudotime)[pseudotime >= 0.8])
day0_lineages <- all_data$assigned_lineage[day0_cells]
tmp <- sort(table(day0_lineages), decreasing = T)
lineage_vec <- names(tmp)
length(lineage_vec)

treatment <- "DABTRAM"
cols.highlight_vec <- c("darkblue", "#DE2D26", "lightpink", "#00CC33")
names(cols.highlight_vec) <- c("week5", "day10_front", "day10_other", "day0")

pdf("../../../../out/figures/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime_frontrunner-losers.pdf", 
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
  cell_names2a <- intersect(cell_names2, day0_cells)
  cell_names2b <- setdiff(cell_names2, day0_cells)
  
  cell_idx3 <- intersect(which(all_data$dataset == "day0"),
                         which(all_data$assigned_lineage == lineage_name))
  cell_names3 <- colnames(all_data)[cell_idx3]
  
  cells.highlight_list <- list(week5 = cell_names,
                               day10_front = cell_names2a,
                               day10_other = cell_names2b,
                               day0 = cell_names3)
  tmp <- cols.highlight_vec
  len <- sapply(cells.highlight_list, length)
  cells.highlight_list <- cells.highlight_list[len > 0]
  tmp <- tmp[len > 0]
  # DOWNBAD BORKED MONKASTEER: The code for Seurat::DimPlot is so scuffed
  # I don't know why the colors are flipping, so I'm doing something super janky to flip it back
  if("day10_other" %in% names(cells.highlight_list)){
    tmp["day10_other"] <- "#DE2D26"
    tmp["day10_front"] <- "lightpink"
  }
  
  p1 <- Seurat::DimPlot(all_data, reduction = "umap",
                        cells.highlight = cells.highlight_list,
                        cols.highlight = tmp)
  p1 <- p1 + ggplot2::ggtitle(paste0("DABTRAM: ", lineage_name, " (", length(cell_names), " week5 cells,\nfrom ",
                                     length(cell_names2a)+length(cell_names2b), " day10 cells (", length(cell_names2a), " in frontrunner position)\n",
                                     "and ", length(cell_idx3), " day0 cells)")) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=9)) + Seurat::NoLegend()
  print(p1)
}

dev.off()