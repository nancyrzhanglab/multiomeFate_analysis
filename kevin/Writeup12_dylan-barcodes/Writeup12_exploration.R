rm(list=ls())
library(Seurat)
library(multiomeFate)
load("~/nzhanglab/data/GoogleDrive_SydneyShafferLab/Dylan_barcode_09-10-2024/all_naive.normalized_named.RData")

all_naive

lin_mat <- SeuratObject::LayerData(all_naive, layer = "counts", assay = "lineage")
dim(lin_mat)

#################

set.seed(10)
res <- multiomeFate:::barcode_clustering(lin_mat,
                                         cell_lower_limit = 15,
                                         cor_threshold = 0.65,
                                         warn_merging = TRUE,
                                         verbose = 1)

quantile(res$minimum_correlation)
sort(res$minimum_correlation, decreasing = FALSE)[1:10]
length(res$lineage_clusters)
table(sapply(res$lineage_clusters, length))

idx <- which(sapply(res$lineage_clusters, length) == 5)
res$minimum_correlation[names(res$lineage_clusters)[idx]]

res2 <- multiomeFate:::barcode_combine(lin_mat,
                                       lineage_clusters = res$lineage_clusters,
                                       verbose = 1)

dim(res2)
quantile(res2@x, probs = seq(0,1,length.out=21))

