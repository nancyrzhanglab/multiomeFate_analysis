rm(list=ls())
# atac_mat <- read.table("../../../../data/Greenleaf_humancortex/GSE162170_multiome_atac_counts.tsv")
# atac_mat <- as.matrix(atac_mat)
# atac_mat <- Matrix::Matrix(atac_mat2, sparse = T)
# 
# atac_peaks <- read.table("../../../../data/Greenleaf_humancortex/GSE162170_multiome_atac_consensus_peaks.txt",
#                          header = T)
# rowname_vec <- apply(atac_peaks, 1, function(x){
#   paste0(x["seqnames"], "_", 
#          gsub(" ", "", as.character(x["start"])), "_", 
#          gsub(" ", "", as.character(x["end"])))
# })
# rownames(atac_mat) <- rowname_vec
# 
# cell_metadata <- read.table("../../../../data/Greenleaf_humancortex/GSE162170_multiome_cell_metadata.txt",
#                             header = T)
# colnames(atac_mat) <- cell_metadata$Cell.ID
# save(atac_mat, file = "../../../../data/Greenleaf_humancortex/GSE162170_multiome_atac_counts.RData")
# 
# rna_mat <- read.table("../../../../data/Greenleaf_humancortex/GSE162170_multiome_rna_counts.tsv")
# rna_mat <- as.matrix(rna_mat)
# rna_mat <- Matrix::Matrix(rna_mat, sparse = T)
# save(rna_mat, file = "../../../../data/Greenleaf_humancortex/GSE162170_multiome_rna_counts.RData")

#####################
# compare/contrast with https://github.com/GreenleafLab/brainchromatin/blob/main/R/Load_Environment.R, but I don't see much from this...
library(Seurat)
load("../../../../data/Greenleaf_humancortex/GSE162170_multiome_atac_counts.RData")
load("../../../../data/Greenleaf_humancortex/GSE162170_multiome_rna_counts.RData")
cell_metadata <- read.table("../../../../data/Greenleaf_humancortex/GSE162170_multiome_cell_metadata.txt",
                            header = T)
cluster_info <- read.table("../../../../data/Greenleaf_humancortex/GSE162170_multiome_cluster_names.txt",
                           header = T, sep = "\t")

# create a new Seurat object
cortex <- Seurat::CreateSeuratObject(counts = rna_mat)
cortex[["ATAC"]] <- Seurat::CreateAssayObject(counts = atac_mat)
celltype <- cell_metadata$seurat_clusters
cluster_info <- cluster_info[which(cluster_info$Assay == "Multiome RNA"),]
for(i in unique(celltype)){
  idx <- which(celltype == i)
  label <- cluster_info[which(cluster_info$Cluster.ID == i), "Cluster.Name"]
  celltype[idx] <- label
}
cortex[["celltype"]] <- celltype

Seurat::DefaultAssay(cortex) <- "RNA"
set.seed(10)
cortex <- Seurat::SCTransform(cortex)
cortex <- Seurat::RunPCA(cortex, verbose = FALSE)
set.seed(10)
cortex <- Seurat::RunUMAP(cortex, dims = 1:50, 
                          reduction.name = 'umap.rna', 
                          reduction.key = 'rnaUMAP_')

Seurat::DefaultAssay(cortex) <- "ATAC"
set.seed(10)
cortex <- Signac::RunTFIDF(cortex)
cortex <- Signac::FindTopFeatures(cortex, min.cutoff = 'q0')
cortex <- Signac::RunSVD(cortex)
set.seed(10)
cortex <- Seurat::RunUMAP(cortex, reduction = 'lsi', dims = 2:50, 
                          reduction.name = "umap.atac", 
                          reduction.key = "atacUMAP_")

set.seed(10)
cortex <- Seurat::FindMultiModalNeighbors(cortex, 
                                          reduction.list = list("pca", "lsi"), 
                                          dims.list = list(1:50, 2:50))
cortex <- Seurat::RunUMAP(cortex, 
                          nn.name = "weighted.nn", 
                          reduction.name = "wnn.umap", 
                          reduction.key = "wnnUMAP_")

save(cortex, file = "../../../../data/Greenleaf_humancortex/09282021_seurat_processed.RData")
