rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../data/Sydney_stressors_2021-09-24/all_data_SCT.RData")

# presumably Lineage is how we track which cells stem from which?

dim(all_data[["RNA"]]@counts)
all_data[["RNA"]]@counts[1:5,1:5]
dim(all_data[["SCT"]]@data)

head(all_data@meta.data)
table(all_data@meta.data$Original_condition)

plot1 <- Seurat::DimPlot(all_data, reduction = "umap", group.by = "Original_condition", 
                         label = TRUE, repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Colored by treatment")
ggplot2::ggsave(filename = "../../../../out/figures/Writeup3d/Writeup3d_sydney_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

zz = all_data[["lineage"]]
zz@counts[1:5,1:5]
mat <- zz@counts
which(mat[1,] != 0)
# idx_list <- lapply(1:nrow(mat), function(j){
#   which(mat[j,] != 0)
# })
mat2 <- Matrix::t(mat)
idx_list <- lapply(1:ncol(mat2), function(j){
  mat2@i[(mat2@p[j]+1):(mat2@p[j+1])]+1
})

len_vec <- sapply(idx_list, length)
quantile(len_vec)
len_vec <- len_vec[len_vec > 2]
num_cell_per_lineage <- len_vec
quantile(num_cell_per_lineage, probs = seq(0,1,length.out=11))

len_vec <- sapply(idx_list, length)
only_one_idx <- which(len_vec == 1)
cell_idx <- unlist(lapply(only_one_idx, function(x){idx_list[[x]]}))
cell_condition <- all_data@meta.data$Original_condition[cell_idx]
which(cell_condition == "acid")[1] # yields 1653
only_one_idx[1653] # yields  2799
which(mat[2799,] != 0) # yields 15319
all_data@meta.data$Original_condition[15319]

idx <- which(len_vec == 435)
cell_idx <- which(mat[idx,] != 0)
table(all_data@meta.data$Original_condition[cell_idx])

mat@x <- rep(1, length(mat@x))
tmp <- sparseMatrixStats::colSums2(mat)
quantile(tmp)
cell_idx <- which(tmp == 12)[1]
which(all_data[["lineage"]]@counts[,201] != 0)
