rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
all_data$common_tcca <- NULL
all_data$distinct1_tcca <- NULL
all_data$distinct2_tcca <- NULL
all_data <- subset(all_data, dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))
n <- ncol(all_data)

Seurat::DefaultAssay(all_data) <- "Saver"
mat_1 <- Matrix::t(all_data[["Saver"]]@data[Seurat::VariableFeatures(object = all_data),])
Seurat::DefaultAssay(all_data) <- "ATAC"
mat_2 <- Matrix::t(all_data[["ATAC"]]@data[Seurat::VariableFeatures(object = all_data),])

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

#################################

set.seed(10)
multiSVD_obj <- tiltedCCA:::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                            dims_1 = 1:50, dims_2 = 2:50,
                                            center_1 = T, center_2 = F,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = T,
                                            rescale_1 = F, rescale_2 = T,
                                            scale_1 = T, scale_2 = F,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = NULL, 
                                           large_clustering_2 = NULL, 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 15,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 15,
                                         verbose = 2)

tmp <- all_data; tmp_mat <- multiSVD_obj$laplacian_list$common_laplacian
colnames(tmp_mat) <- paste0("tmp_", 1:ncol(tmp_mat))
set.seed(10); tmp_umap <- Seurat::RunUMAP(tmp_mat)@cell.embeddings
tmp_umap_full <- matrix(NA, nrow = ncol(tmp), ncol = 2)
for(i in 1:length(multiSVD_obj$metacell_obj$metacell_clustering_list)){
  idx <- multiSVD_obj$metacell_obj$metacell_clustering_list[[i]]
  tmp_umap_full[idx,] <- rep(tmp_umap[i,], each = length(idx))
}
set.seed(10)
tmp_umap_full <- jitter(tmp_umap_full)
rownames(tmp_umap_full) <- colnames(tmp)
tmp[["common_laplacian"]] <- Seurat::CreateDimReducObject(tmp_umap_full, key = "commonLapUMAP")
plot1 <- Seurat::DimPlot(tmp, reduction = "common_laplacian",
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Target common laplacian"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6c/Writeup6c_tcca_RNA-ATAC_COCL2_laplacian.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##################3

multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = F)

save(multiSVD_obj, all_data,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_DABTRAM.RData")

#################################

set.seed(10)
all_data[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                          what = "common",
                                                          aligned_umap_assay = "umap",
                                                          seurat_obj = all_data,
                                                          seurat_assay = "RNA",
                                                          verbose = 1)
set.seed(10)
all_data[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                             what = "distinct_1",
                                                             aligned_umap_assay = "umap",
                                                             seurat_obj = all_data,
                                                             seurat_assay = "RNA",
                                                             verbose = 1)
set.seed(10)
all_data[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                             what = "distinct_2",
                                                             aligned_umap_assay = "umap",
                                                             seurat_obj = all_data,
                                                             seurat_assay = "RNA",
                                                             verbose = 1)

save(multiSVD_obj, all_data,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6c/Writeup6c_tcca_RNA-ATAC_COCL2.RData")
