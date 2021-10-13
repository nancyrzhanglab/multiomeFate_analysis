rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/kevin/Writeup3d/09302021_sydney_de_preprocess.RData")
ls()

lin_mat <- all_data[["lineage"]]@counts
tmp <- Matrix::t(lin_mat)
lin_idx_list <- lapply(1:ncol(tmp), function(j){
  .nonzero_col(tmp, j)
})
names(lin_idx_list) <- colnames(tmp)
factor_vec <- as.factor(all_data@meta.data$Original_condition)
tabulate_mat <- t(sapply(lin_idx_list, function(idx){
  table(factor_vec[idx])
}))
rownames(tabulate_mat) <- colnames(tmp)
tabulate_mat <- tabulate_mat[order(tabulate_mat[,"naive"], decreasing = T),]

############

# first plot, all the cells in the dataset
plot1 <- Seurat::DimPlot(all_data, reduction = "umap", group.by = "Original_condition", 
                         label = TRUE, repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Colored by treatment")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3d/Writeup3d_sydney_de_all.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


# next, focusing on Cocl3-survival, plot the naive cells
celltypes <- unique(all_data@meta.data$Original_condition)
cell_terminal <- which(all_data@meta.data$Original_condition == "cocl2")
lineage_terminal <- which(sparseMatrixStats::rowSums2(lin_mat[,cell_terminal,drop=F]) > 0)
tmp <- Matrix::t(lin_mat)
all_cell_names <- unique(unlist(lapply(lineage_terminal, function(j){
  colnames(lin_mat)[.nonzero_col(tmp, j)]
})))
only_naive <- rownames(all_data@meta.data)[which(all_data@meta.data$Original_condition == "naive")]
naive_terminal <- all_cell_names[which(all_cell_names %in% only_naive)]
length(naive_terminal)

n <- ncol(all_data[["SCT"]])
naive_bool <- rep(NA, n)
names(naive_bool) <- colnames(all_data)
naive_bool[only_naive] <- 0
naive_bool[naive_terminal] <- 1
table(naive_bool)
xlim <- range(all_data[["umap"]]@cell.embeddings[,1])
ylim <- range(all_data[["umap"]]@cell.embeddings[,2])

all_data[["naive"]] <- naive_bool
cells <- unique(c(only_naive, naive_terminal))
plot1 <- Seurat::DimPlot(all_data, reduction = "umap", group.by = "naive", 
                         label = TRUE, repel = TRUE, label.size = 2.5,
                         cells = cells) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Colored by Cocl3-survival")
plot1 <- plot1 + ggplot2::expand_limits(x = xlim, y = ylim)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3d/Writeup3d_sydney_de_cocl3_survival.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

#########################

# some lineages to try
# Lin611467 (multiple fates), Lin283219 (dies to all), Lin156893 and Lin206249 (just cocl3 and naive)

col_palette <- scales::hue_pal()(6)
lineage_vec <- c("Lin611467", "Lin283219", "Lin156893", "Lin206249")
for(lineage in lineage_vec){
  cells <- which(lin_mat[lineage,] != 0)
  col_vec <- col_palette[which(sort(unique(all_data@meta.data$Original_condition)) %in% all_data@meta.data$Original_condition[cells])]
  plot1 <- Seurat::DimPlot(all_data, reduction = "umap", group.by = "Original_condition", 
                           label = TRUE, repel = TRUE, label.size = 2.5,
                           cells = cells, cols = col_vec) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Lineage ", lineage, "\nColored by treatment"))
  plot1 <- plot1 + ggplot2::expand_limits(x = xlim, y = ylim)
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_sydney_de_", lineage, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}


