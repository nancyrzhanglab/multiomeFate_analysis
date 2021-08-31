rm(list=ls())

library(Seurat); library(Signac)
source("diffusion_functions.R")
source("metacell_construction.R")
source("metacell_graph_plot.R")

load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed_de.RData")

Seurat::DefaultAssay(mbrain3) <- "ATAC"
mbrain3 <- Seurat::FindNeighbors(
  object = mbrain3,
  reduction = 'lsi',
  dims = 2:50
)

set.seed(10)
mbrain3 <- Seurat::FindClusters(
  object = mbrain3,
  algorithm = 3,
  resolution = 10,
  verbose = T
)

#######################

clustering <- as.character(mbrain3@meta.data$ATAC_snn_res.10)
metacell_mat <- form_metacell_matrix(mbrain3[["lsi"]]@cell.embeddings, clustering)
adj_mat <- form_snn_graph(metacell_mat, k = 10, distance = "cosine")
P <- form_transition(adj_mat, normalize = T)
res <- extract_eigen(P, check = T)

target_vec <- c("4", # radial
                "46", # oligo
                "11", # forebrain
                "1", "8", #cortical1, cortical2
                "31", "5", # neuroblast1, neuroblast2
                "27") #glio
target_name_vec <- c("Radial", "Oligo", "Forebrain", 
                     "Cortical1", "Cortical2", "Neuro1", "Neuro2",
                     "Glio")

uniq_clust <- sort(unique(clustering))
dist_mat <- sapply(target_vec, function(target){
  i <- which(rownames(metacell_mat) == target)
  sapply(1:nrow(metacell_mat), function(j){
    diffusion_distance(res$eigenvalues, res$right_vector,
                       i, j, time_vec = c(1:3))
  })
})
rownames(dist_mat) <- rownames(metacell_mat)

png(paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_metacell_knn_diffusion1.png"),
    height = 4500, width = 4500, res = 300, units = "px")
par(mfrow = c(3,3))
for(i in 1:length(target_vec)){
  plot_metacell_graph(mbrain3[["wnn.umap"]]@cell.embeddings,
                      clustering,
                      adj_mat,
                      feature_vec = dist_mat[,i],
                      zlim = range(dist_mat[dist_mat > 1e-6]),
                      xlab = "wnnUMAP_1",
                      ylab = "wnnUMAP_2",
                      main = paste0("(Norm.) Diffusion dist to ", target_name_vec[i]))
}
plot_legend(zlim = range(dist_mat[dist_mat > 1e-6]),
            offset = -0.75)
graphics.off()

###############################

tmp <- rbind(1:3, function(i){
  
})

