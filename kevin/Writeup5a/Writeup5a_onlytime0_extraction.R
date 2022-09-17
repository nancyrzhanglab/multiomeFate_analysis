rm(list=ls())
library(Seurat)
library(Signac)

load("../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging_complete.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

time0 <- all_data
time0[["geneActivity"]] <- NULL
time0[["Lineage"]] <- NULL
time0[["spliced"]] <- NULL
time0[["unspliced"]] <- NULL
time0[["activity.umap"]] <- NULL
time0[["fasttopic_CIS"]] <- NULL
time0[["fasttopic_COCL2"]] <- NULL
time0[["fasttopic_DABTRAM"]] <- NULL
time0[["pca"]] <- NULL
time0[["lsi"]] <- NULL
time0[["umap"]] <- NULL
time0[["atac.umap"]] <- NULL
time0[["saverpca"]] <- NULL
time0[["saverumap"]] <- NULL

keep_vec <- rep(0, ncol(time0))
keep_vec[which(time0$dataset == "day0")] <- 1
time0$keep <- keep_vec
time0 <- subset(time0, keep == 1)

time0$nCount_geneActivity <- NULL
time0$nFeature_geneActivity <- NULL
time0$nCount_Saver <- NULL
time0$nFeature_Saver <- NULL
time0$nCount_Lineage <- NULL
time0$nFeature_Lineage <- NULL
time0$nCount_spliced <- NULL
time0$nFeature_spliced <- NULL
time0$nCount_unspliced <- NULL
time0$nFeature_unspliced <- NULL

save(time0, date_of_run, session_info,
     file = "../../../out/kevin/Writeup5a/Writeup5a_time0Only.RData")




