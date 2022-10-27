rm(list=ls())

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tiltedCCA)

load("../../../../out/kevin/Writeup6/Writeup6_tcca_selected-genes.RData")
load("../../../../out/kevin/Writeup5a/Writeup5a_tcca_RNA-geneActivity.RData")
all_data2 <- all_data
load("../../../../out/kevin/Writeup6/Writeup6_all-data_lineage-assigned.RData")
source("../Writeup5a/color_palette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#####################

multiSVD_obj[["common_mat_1"]] <- NULL
multiSVD_obj[["distinct_mat_1"]] <- NULL
multiSVD_obj[["common_dimred_2"]] <- NULL
multiSVD_obj[["distinct_dimred_2"]] <- NULL

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

#####################

common_mat_1 <- multiSVD_obj$common_mat_1[,selection_res$selected_variables]
common_mat_2 <- multiSVD_obj$common_mat_2[,selection_res$selected_variables]

common_mat_1 <- tiltedCCA:::.normalize_svd(input_obj = common_mat_1,
                                           averaging_mat = NULL,
                                           center = T,
                                           dims = 1:50,
                                           normalize_row = T,
                                           normalize_singular_value = T,
                                           recenter = F,
                                           rescale = F,
                                           scale = T,
                                           scale_max = NULL)
common_mat_2 <- tiltedCCA:::.normalize_svd(input_obj = common_mat_2,
                                           averaging_mat = NULL,
                                           center = T,
                                           dims = 2:50,
                                           normalize_row = T,
                                           normalize_singular_value = T,
                                           recenter = F,
                                           rescale = F,
                                           scale = T,
                                           scale_max = NULL)
common_mat <- cbind(common_mat_1, common_mat_2)

set.seed(10)
seurat_umap <- Seurat::RunUMAP(common_mat, 
                               assay = "RNA",
                               verbose = T)
rownames(seurat_umap@cell.embeddings) <- rownames(common_mat)
all_data[["common_tcca"]] <- seurat_umap

umap_mat <- all_data[["common_tcca"]]@cell.embeddings
dataset_vec <- all_data$dataset
set.seed(10)
idx <- sample(1:nrow(umap_mat))
umap_mat <- umap_mat[idx,]
dataset_vec <- dataset_vec[idx]
count_mat <- all_data[["RNA"]]@counts[,idx]

all_data3 <- Seurat::CreateSeuratObject(counts = count_mat)
all_data3$dataset <- dataset_vec
all_data3[["common_tcca"]] <- Seurat::CreateDimReducObject(embeddings = umap_mat)

plot1 <- Seurat::DimPlot(all_data3, reduction = "common_tcca",
                         group.by = "dataset", 
                         cols = col_palette, pt.size = .1)
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6/Writeup6_tcca-common_RNA-geneActivity_100selected-genes.png"),
                plot1, device = "png", width = 7, height = 6, units = "in")

#################

umap_mat <- all_data2[["common_tcca"]]@cell.embeddings
dataset_vec <- all_data2$dataset
set.seed(10)
idx <- sample(1:nrow(umap_mat))
umap_mat <- umap_mat[idx,]
dataset_vec <- dataset_vec[idx]
count_mat <- all_data2[["RNA"]]@counts[,idx]

all_data3 <- Seurat::CreateSeuratObject(counts = count_mat)
all_data3$dataset <- dataset_vec
all_data3[["common_tcca"]] <- Seurat::CreateDimReducObject(embeddings = umap_mat)

plot1 <- Seurat::DimPlot(all_data3, reduction = "common_tcca",
                         group.by = "dataset", 
                         cols = col_palette, pt.size = .1)
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6/Writeup6_tcca-common_RNA-geneActivity_full-genes.png"),
                plot1, device = "png", width = 7, height = 6, units = "in")

##################

umap_mat <- all_data[["saverumap"]]@cell.embeddings
dataset_vec <- all_data$dataset
set.seed(10)
idx <- sample(1:nrow(umap_mat))
umap_mat <- umap_mat[idx,]
dataset_vec <- dataset_vec[idx]
count_mat <- all_data[["RNA"]]@counts[,idx]

all_data3 <- Seurat::CreateSeuratObject(counts = count_mat)
all_data3$dataset <- dataset_vec
all_data3[["saverumap"]] <- Seurat::CreateDimReducObject(embeddings = umap_mat)

plot1 <- Seurat::DimPlot(all_data3, reduction = "saverumap",
                         group.by = "dataset", 
                         cols = col_palette, pt.size = .1)
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Saver UMAP 2") + ggplot2::xlab("Saver UMAP 1")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6/Writeup6_Saver-umap_full-genes.png"),
                plot1, device = "png", width = 7, height = 6, units = "in")


##################

umap_mat <- all_data[["activity.umap"]]@cell.embeddings
dataset_vec <- all_data$dataset
set.seed(10)
idx <- sample(1:nrow(umap_mat))
umap_mat <- umap_mat[idx,]
dataset_vec <- dataset_vec[idx]
count_mat <- all_data[["RNA"]]@counts[,idx]

all_data3 <- Seurat::CreateSeuratObject(counts = count_mat)
all_data3$dataset <- dataset_vec
all_data3[["activity.umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_mat)

plot1 <- Seurat::DimPlot(all_data3, reduction = "activity.umap",
                         group.by = "dataset", 
                         cols = col_palette, pt.size = .1)
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Saver UMAP 2") + ggplot2::xlab("Saver UMAP 1")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6/Writeup6_geneActivity-umap_full-genes.png"),
                plot1, device = "png", width = 7, height = 6, units = "in")

#########################

count_mat <- all_data[["Saver"]]@scale.data[selection_res$selected_variables,]
set.seed(10)
svd_res <- irlba::irlba(count_mat, nv = 50)
dimred <- tiltedCCA:::.mult_mat_vec(svd_res$v, svd_res$d)

set.seed(10)
seurat_umap <- Seurat::RunUMAP(dimred, 
                               assay = "RNA",
                               verbose = T)
rownames(seurat_umap@cell.embeddings) <- colnames(all_data)
all_data[["Sumap"]] <- seurat_umap

umap_mat <- all_data[["Sumap"]]@cell.embeddings
dataset_vec <- all_data$dataset
set.seed(10)
idx <- sample(1:nrow(umap_mat))
umap_mat <- umap_mat[idx,]
dataset_vec <- dataset_vec[idx]
count_mat <- all_data[["RNA"]]@counts[,idx]

all_data3 <- Seurat::CreateSeuratObject(counts = count_mat)
all_data3$dataset <- dataset_vec
all_data3[["Sumap"]] <- Seurat::CreateDimReducObject(embeddings = umap_mat)

plot1 <- Seurat::DimPlot(all_data3, reduction = "Sumap",
                         group.by = "dataset", 
                         cols = col_palette, pt.size = .1)
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Saver UMAP 2") + ggplot2::xlab("Saver UMAP 1")
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6/Writeup6_Saver-umap_100selected-genes.png"),
                plot1, device = "png", width = 7, height = 6, units = "in")

################################

# now let's try to plot the cells in a lineage
rm(list=c("all_data2", "all_data3", "common_mat", "common_mat_1", "common_mat_2", 
          "count_mat", "dataset_vec", "dimred", "seurat_umap", "svd_res", "umap_mat"))
gc(T)

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

plotting_func <- function(lineage_name,
                          seurat_obj,
                          dim_name){
  par(mfrow = c(3,3), mar = c(0.5, 0.5, 4, 0.5))
  stopifnot(lineage_name %in% seurat_obj$assigned_lineage)
  
  treatment_vec <- c("day0", "day10_CIS", "week5_CIS", 
                     NA, "day10_COCL2", "week5_COCL2",
                     NA, "day10_DABTRAM", "week5_DABTRAM")
  for(kk in treatment_vec){
    if(is.na(kk)){
      graphics::plot(NULL, type = "n",  xlim = c(0,1), ylim = c(0, 1),
                     xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", 
                     xlab = '', ylab = '', main = "", frame=FALSE)
    } else {
      umap_mat <- all_data[[dim_name]]@cell.embeddings
      idx <- which(all_data$dataset == kk)
      graphics::plot(umap_mat[idx,1], umap_mat[idx,2], 
                     pch = 16, cex = 1, col = "gray",
                     xlab = "", ylab = "", main = kk,
                     xaxt = "n", yaxt = "n")
      idx <- which(all_data$dataset == kk & all_data$assigned_lineage == lineage_name)
      graphics::points(umap_mat[idx,1], umap_mat[idx,2],
                       pch = 16, cex = 1.8, col = "white")
      graphics::points(umap_mat[idx,1], umap_mat[idx,2],
                       pch = 16, cex = 1.5, col = rgb(1,0,0,0.2))
    }
  }
}


lineage_vec <- c("Lin10837", "Lin33594", "Lin43542", #CIS
                 "Lin51840", "Lin1550", # COCL2
                 "Lin52442", "Lin112532", "Lin33330", "Lin1329") #DABTRAM
for(lineage_name in lineage_vec){
  print(lineage_name)
  png(paste0("../../../../out/figures/Writeup6/Writeup6_tcca-common_RNA-geneActivity_100selected-genes_lineage-", lineage_name, ".png"),
      height = 3000, width = 3000, units = "px", res = 300)
  plotting_func(lineage_name = lineage_name,
                seurat_obj = all_data,
                dim_name = "common_tcca")
  graphics.off()
}

for(lineage_name in lineage_vec){
  print(lineage_name)
  png(paste0("../../../../out/figures/Writeup6/Writeup6_Saver-umap_100selected-genes_lineage-", lineage_name, ".png"),
      height = 3000, width = 3000, units = "px", res = 300)
  plotting_func(lineage_name = lineage_name,
                seurat_obj = all_data,
                dim_name = "Sumap")
  graphics.off()
}
