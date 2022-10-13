rm(list=ls())
library(Seurat); library(Signac)

load("../../../../out/kevin/Writeup4e/Writeup4e_tcca_RNA-ATAC.RData")
source("color_palette.R")

umap_mat <- all_data[["common_tcca"]]@cell.embeddings
dataset_vec <- all_data$dataset
set.seed(10)
idx <- sample(1:nrow(umap_mat))
umap_mat <- umap_mat[idx,]
dataset_vec <- dataset_vec[idx]
count_mat <- all_data[["RNA"]]@counts[,idx]

all_data2 <- Seurat::CreateSeuratObject(counts = count_mat)
all_data2$dataset <- dataset_vec
all_data2[["common_tcca"]] <- Seurat::CreateDimReducObject(embeddings = umap_mat)

plot1 <- Seurat::DimPlot(all_data2, reduction = "common_tcca",
                         group.by = "dataset", 
                         cols = col_palette, pt.size = .3)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup5a/Writeup5a_tcca_RNA-ATAC_common-slides.png"),
                plot1, device = "png", width = 4.5, height = 2.5, units = "in")


