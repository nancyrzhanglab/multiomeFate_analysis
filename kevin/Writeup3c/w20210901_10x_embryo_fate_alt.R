rm(list=ls())

load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed_de.RData")

library(Seurat); library(Signac); library(multiomeFate)

load_func <- function(){
  source("candidate_method_alt.R")
  source("chromatin_potential_alt.R")
  source("chromatin_potential_preparation_alt.R")
  source("diffusion_functions.R")
  source("estimation_method_alt.R")
  source("forming_method_alt.R")
  source("matching_method.R")
  source("metacell_construction.R")
  source("metacell_graph_plot.R")
}
load_func()

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

###############################

clustering <- as.character(mbrain3@meta.data$ATAC_snn_res.10)
metacell_mat <- form_metacell_matrix(mbrain3[["lsi"]]@cell.embeddings, clustering)
snn <- form_snn_graph(metacell_mat, k = 10, distance = "cosine")
P <- form_transition(snn, normalize = T,
                     lazy_param = 0.85,
                     teleport_param = 0.99)
res <- extract_eigen(P, check = T)
n <- nrow(metacell_mat)
diffusion_dist <- sapply(1:n, function(i){
  sapply(1:n, function(j){
    diffusion_distance(res$eigenvalues, res$right_vector,
                       i, j)
  })
})
rownames(diffusion_dist) <- rownames(metacell_mat)
colnames(diffusion_dist) <- rownames(metacell_mat)

###############################

# grab the relevant genes and peaks
mat_x <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))
peak_names <- sort(unique(mbrain3[["ATAC"]]@links$peak))
gene_names <- sort(unique(mbrain3[["ATAC"]]@links$gene))
mat_x <- mat_x[,which(colnames(mat_x) %in% peak_names)]
mat_y <- mat_y[,which(colnames(mat_y) %in% gene_names)]
mat_x <- form_metacell_matrix(mat_x, clustering, func = mean)
mat_y <- form_metacell_matrix(mat_y, clustering, func = mean)

# create the hash map
p1 <- ncol(mat_x); p2 <- ncol(mat_y); n <- nrow(mat_x)
ht_map <- hash::hash()
for(i in 1:p2){
  gene_name <- colnames(mat_y)[i]
  idx <- which(mbrain3[["ATAC"]]@links$gene == gene_name)
  peak_names <- mbrain3[["ATAC"]]@links$peak[idx]
  ht_map[[as.character(i)]] <- which(colnames(mat_x) %in% peak_names)
}

vec_start <- c("4", "48")
list_end <- list(c("28", "36", "19", "11"), #forebrain
                 c("46"), # oligo
                 c("52", "6", "0", "8"), #cortical2
                 c("1", "2", "3", "7", #cortical1
                   "14", "18", 
                   "24", "25", "29", 
                   "30", "32", "37", "39",
                   "40", "42","44", "47",
                   "53"))

# check the assignments
length(intersect(which(mbrain3@meta.data$ATAC_snn_res.10 %in% vec_start),
                 which(mbrain3@meta.data$new_seurat_clusters == 15)))/
  length(which(mbrain3@meta.data$new_seurat_clusters == 15))
terminal_list <- list("6", "16", "9", c("1", "2", "4"))
sapply(1:length(list_end), function(i){
  length(intersect(which(mbrain3@meta.data$ATAC_snn_res.10 %in% list_end[[i]]),
            which(mbrain3@meta.data$new_seurat_clusters %in% terminal_list[[i]])))/
  length(which(mbrain3@meta.data$new_seurat_clusters %in% terminal_list[[i]]))
})

rm(list = c("gene_name", "gene_names", "i", "idx",
            "metacell_mat", "n", "P",
            "p1", "p2", "peak_names", "res", "terminal_list"))
rm(list = ls()[!ls() %in% c("diffusion_dist", 
                            "mat_x", "mat_y",
                            "snn", "vec_start",
                            "list_end", "ht_map", "load_func",
                            "de_combined")])
load_func()

#############################

prep_obj <- chromatin_potential_prepare2(mat_x, 
                                         mat_y, 
                                         snn,
                                         diffusion_dist,
                                         vec_start, 
                                         list_end,
                                         ht_map)


