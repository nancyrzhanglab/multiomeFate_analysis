rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
source("snn_helpers.R")

load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_processed2.RData")


set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

rna_dimred <- all_data[["pca"]]@cell.embeddings[,1:30]
l2_vec <- apply(rna_dimred, 1, function(x){sqrt(sum(x^2))})
l2_vec <- pmax(l2_vec, 1e-3)
rna_dimred <- .mult_vec_mat(1/l2_vec, rna_dimred)
rna_snn <- .form_snn_mat(rna_dimred, 
                        num_neigh = 30,
                        bool_cosine = T, 
                        bool_intersect = F)

atac_dimred <- all_data[["lsi"]]@cell.embeddings[,2:30]
l2_vec <- apply(atac_dimred, 1, function(x){sqrt(sum(x^2))})
l2_vec <- pmax(l2_vec, 1e-3)
atac_dimred <- .mult_vec_mat(1/l2_vec, atac_dimred)
atac_nn <- .form_snn_mat(atac_dimred, 
                         num_neigh = 30,
                         bool_cosine = T, 
                         bool_intersect = F)

n <- nrow(rna_dimred)
jaccard_vec <- sapply(1:n, function(i){
  vec1 <- rna_nn$nn.idx[i,-1]
  vec2 <- atac_nn$nn.idx[i,-1]
  
  length(intersect(vec1, vec2))/length(unique(c(vec1, vec2)))
})
quantile(jaccard_vec)
all_data$jaccard <- jaccard_vec

plot1 <-Seurat::FeaturePlot(all_data, reduction = "umap",
                            features = "jaccard")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Jaccard similarity\n(k=30), RNA UMAP"))

plot2 <-Seurat::FeaturePlot(all_data, reduction = "atac.umap",
                            features = "jaccard")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Jaccard similarity\n(k=30), ATAC UMAP"))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_jaccard_umap.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")
