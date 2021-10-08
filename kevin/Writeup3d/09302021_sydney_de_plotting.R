rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/kevin/Writeup3d/09302021_sydney_de.RData")
ls()

#################

plot1 <- Seurat::DimPlot(all_data, reduction = "umap", group.by = "Original_condition", 
                         label = TRUE, repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Colored by treatment")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3d/Writeup3d_sydney_de_all.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

n <- ncol(all_data[["SCT"]])
celltypes <- unique(all_data@meta.data$Original_condition)
terminal_celltypes <- setdiff(celltypes, "naive")

naive_lineage <- rep(NA, n)
for(celltype in terminal_celltypes){
  cells <- which(all_data@meta.data$Original_condition == celltype)
  lineage_idx <- which(sparseMatrixStats::rowSums2(all_data[["lineage"]]@counts[,cells]) > 0)
  lineages <- rownames(all_data[["lineage"]]@counts)[lineage_idx]
  stopifnot(all(lineages %in% keep_lin_list[[which(names(keep_lin_list) == celltype)]]))
  
  lineage_cell_idx1 <- which(sparseMatrixStats::colSums2(all_data[["lineage"]]@counts[lineages,]) > 0)
  lineage_cell_idx2 <- which(all_data@meta.data$Original_condition == "naive")
  idx <- intersect(lineage_cell_idx1, lineage_cell_idx2)
  naive_lineage[idx] <- celltype
}
all_data[["naive"]] <- naive_lineage
col_vec <- scales::hue_pal()(length(celltypes))
col_vec <- col_vec[-5]

plot1 <- Seurat::DimPlot(all_data, reduction = "umap", group.by = "naive", 
                         label = TRUE, repel = TRUE, label.size = 2.5, cols = col_vec) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Colored by naive fate")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3d/Writeup3d_sydney_de_naive.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


