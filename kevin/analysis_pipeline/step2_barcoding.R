rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

load("../../../../out/kevin/analysis_pipeline/step1_peakmerging.RData")

lin_mat <- all_data[["Lineage"]]@counts
lin_mat <- lin_mat[Matrix::rowSums(lin_mat) > 0,]

res_clusters <- multiomeFate::barcode_clustering(lin_mat)
cluster_size <- sapply(res_clusters$lineage_clusters, length)

# combining lineages
lin_mat <- as.matrix(lin_mat)
lin_mat2 <- lin_mat
lin_mat <- barcode_combine(lin_mat = lin_mat,
                           lineage_clusters = res_clusters$lineage_clusters,
                           verbose = 1)

#########################

posterior_res <- multiomeFate::barcoding_posterior(lin_mat = lin_mat,
                                                   verbose = 1)

maxBhat <- apply(posterior_res$posterior_mat, 2, max) #chosen lineage 
assignment_vec <- multiomeFate::barcoding_assignment(posterior_mat = posterior_res$posterior_mat,
                                                     difference_val = 0.2,
                                                     verbose = 1)

all_data[["Lineage"]]@data <- posterior_res$posterior_mat
all_data$assigned_lineage <- assignment_vec
all_data$assigned_posterior <- maxBhat

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/analysis_pipeline/step2_barcoding.RData")

