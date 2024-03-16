rm(list=ls())
library(Seurat)
library(scCustomize)

##### Create seurat obj #####

# download files from https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation

# load in the metadata
cell_metadata_df <- read.csv("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_metadata.txt.gz",
                             sep = "\t")
rowname_vec <- sapply(1:nrow(cell_metadata_df), function(i){
  paste0(cell_metadata_df[i,"Cell.barcode"], "_", i)
})
rownames(cell_metadata_df) <- rowname_vec

pseudotime_df <- read.csv("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_neutrophil_pseudotime.txt.gz",
                          sep = "\t")
pseudotime_vec <- rep(NA, length(rowname_vec))
pseudotime_vec[1+pseudotime_df$Cell.index] <- pseudotime_df$pseudotime
trajectory_df <- read.csv("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_neutrophil_monocyte_trajectory.txt.gz",
                          sep = "\t")
trajectory_vec <- rep(FALSE, length(rowname_vec))
trajectory_vec[trajectory_df$Cell.index+1] <- TRUE

cell_metadata_df$pseudotime <- pseudotime_vec
cell_metadata_df$trajectory <- trajectory_vec

# create the count matrix
expression_matrix <- Seurat::ReadMtx(mtx = "~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_normed_counts.mtx.gz",
                                     cells = "~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_metadata.txt.gz",
                                     features = "~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_gene_names.txt.gz",
                                     cell.column = 2,
                                     feature.column = 1,
                                     mtx.transpose = T,
                                     skip.cell = 1)
colnames(expression_matrix) <- rowname_vec

# create the barcode matrix
barcode_mat <- Matrix::readMM("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_clone_matrix.mtx.gz")
barcode_mat <- Matrix::t(barcode_mat)
colnames(barcode_mat) <- rowname_vec
rownames(barcode_mat) <- paste0("Lineage_", 1:nrow(barcode_mat))

#########################
seurat_object <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                            data = expression_matrix,
                                            meta.data = cell_metadata_df)
seurat_object[["Lineage"]] <- Seurat::CreateAssayObject(counts = barcode_mat)

Seurat::DefaultAssay(seurat_object) <- "RNA"
set.seed(10)
seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                              selection.method = "vst",
                                              nfeatures = 2000)
seurat_object <- Seurat::ScaleData(seurat_object)
seurat_object <- Seurat::RunPCA(seurat_object,
                                features = Seurat::VariableFeatures(object = seurat_object),
                                verbose = F)
seurat_object <- Seurat::RunUMAP(seurat_object,
                                 dims = 1:30)

# also include the author-generated SPRING-embedding
spring_embedding <- cbind(spring_1 = seurat_object$SPRING.x,
                          spring_2 = seurat_object$SPRING.y)
seurat_object[["SPRING"]] <- Seurat::CreateDimReducObject(embeddings = spring_embedding)

# plot based on cell types
Seurat::Idents(seurat_object) <- "Cell.type.annotation"
p1 <- Seurat::DimPlot(seurat_object,
                      reduction = "umap")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
p2 <- Seurat::DimPlot(seurat_object,
                      reduction = "SPRING")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
plot_all <- cowplot::plot_grid(p1, p2)
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_larry_umap-celltype.png"),
                plot_all, device = "png", width = 10, height = 5, units = "in")


# plot based on time points
Seurat::Idents(seurat_object) <- "Time.point"
p1 <- Seurat::DimPlot(seurat_object,
                      reduction = "umap")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
p2 <- Seurat::DimPlot(seurat_object,
                      reduction = "SPRING")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
plot_all <- cowplot::plot_grid(p1, p2)
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_larry_umap-timepoint.png"),
                plot_all, device = "png", width = 10, height = 5, units = "in")


# plot based on provided trajectory
p1 <- Seurat::DimPlot(seurat_object,
                      group.by = "trajectory",
                      reduction = "umap")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
p2 <- Seurat::DimPlot(seurat_object,
                      group.by = "trajectory",
                      reduction = "SPRING")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
plot_all <- cowplot::plot_grid(p1, p2)
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_larry_umap-trajectory.png"),
                plot_all, device = "png", width = 10, height = 5, units = "in")


# plot based on provided pseudotime
p1 <- scCustomize::FeaturePlot_scCustom(seurat_object,
                                        features = "pseudotime",
                                        reduction = "umap")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
p2 <- scCustomize::FeaturePlot_scCustom(seurat_object,
                                        features = "pseudotime",
                                        reduction = "SPRING")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 8)) 
plot_all <- cowplot::plot_grid(p1, p2)
ggplot2::ggsave(filename = paste0("~/project/Multiome_fate/out/figures/Writeup8/Writeup8_larry_umap-pseudotime.png"),
                plot_all, device = "png", width = 10, height = 5, units = "in")

###########

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Basic formation of Seurat object from the LARRY dataset in https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation")

save(seurat_object,
     date_of_run, session_info, note,
     file = "~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset.RData")

print("Done! :)")

