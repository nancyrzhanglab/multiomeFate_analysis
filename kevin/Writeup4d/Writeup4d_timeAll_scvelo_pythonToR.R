rm(list=ls())
library(hdf5r)
library(loomR)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

lfile <- loomR::connect(filename = "/home/stat/kevinl1/project/Multiome_fate/out/kevin/Writeup4c/2022_04_week5_CIS_scvelo.loom", 
                        mode = "r+", skip.validate = TRUE)
lfile
lfile[["col_attrs"]]
umap_mat <- t(lfile[["col_attrs/X_umap"]][1:2,])
cellID <- lfile[["col_attrs/CellID"]][]
geneName <- lfile[["row_attrs/Gene"]][]
clusters <- lfile[["col_attrs/Clusters"]][]

save(umap_mat, cellID, geneName, clusters, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_scvelo_umap_CIS.RData")

lfile[["layers"]]
spliced_mat <- lfile[["layers/spliced"]][,]
save(spliced_mat, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_spliced_mat_CIS.RData")

unspliced_mat <- lfile[["layers/unspliced"]][,]
save(unspliced_mat, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_unspliced_mat_CIS.RData")
lfile$close_all()

###############################

lfile <- loomR::connect(filename = "/home/stat/kevinl1/project/Multiome_fate/out/kevin/Writeup4c/2022_04_week5_COCL2_scvelo.loom", 
                        mode = "r+", skip.validate = TRUE)

umap_mat <- t(lfile[["col_attrs/X_umap"]][1:2,])
cellID <- lfile[["col_attrs/CellID"]][]
geneName <- lfile[["row_attrs/Gene"]][]
clusters <- lfile[["col_attrs/Clusters"]][]
save(umap_mat, cellID, geneName, clusters, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_scvelo_umap_COCL2.RData")

spliced_mat <- lfile[["layers/spliced"]][,]
save(spliced_mat, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_spliced_mat_COCL2.RData")

unspliced_mat <- lfile[["layers/unspliced"]][,]
save(unspliced_mat, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_unspliced_mat_COCL2.RData")
lfile$close_all()

###############################

lfile <- loomR::connect(filename = "/home/stat/kevinl1/project/Multiome_fate/out/kevin/Writeup4c/2022_04_week5_DABTRAM_scvelo.loom", 
                        mode = "r+", skip.validate = TRUE)

umap_mat <- t(lfile[["col_attrs/X_umap"]][1:2,])
cellID <- lfile[["col_attrs/CellID"]][]
geneName <- lfile[["row_attrs/Gene"]][]
clusters <- lfile[["col_attrs/Clusters"]][]
save(umap_mat, cellID, geneName, clusters, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_scvelo_umap_DABTRAM.RData")

spliced_mat <- lfile[["layers/spliced"]][,]
save(spliced_mat, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_spliced_mat_DABTRAM.RData")

unspliced_mat <- lfile[["layers/unspliced"]][,]
save(unspliced_mat, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_unspliced_mat_DABTRAM.RData")
lfile$close_all()