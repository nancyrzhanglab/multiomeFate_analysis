rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_peakmerging.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_scvelo_umap_COCL2.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_spliced_mat_COCL2.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_unspliced_mat_COCL2.RData")
maestro_day0 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time0_202204/MAESTRO_arc0_rpmatrix.RDS")
maestro_day10 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time10_cocl2_202204/MAESTRO_arc10_cocl2_rpmatrix.RDS")
maestro_week5 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/week5_cocl2_202204/MAESTRO_week5_cocl2_rpmatrix.RDS")

library(Seurat)
library(Signac)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_COCL2", "week5_COCL2"))] <- 1
all_data$keep <- keep_vec
all_data_COCL2 <- subset(all_data, keep == 1)

#####################

# format the spliced and unspliced mat in
spliced_mat <- Matrix::Matrix(t(spliced_mat), sparse = T)
unspliced_mat <- Matrix::Matrix(t(unspliced_mat), sparse = T)
cellID2 <- sapply(cellID, function(str){
  tmp <- strsplit(str, split = ":")[[1]]
  tmp2 <- strsplit(tmp, split = "_")[[1]]
  start <- tmp2[3]
  if(start != "day0"){
    start <- paste0(start, "_COCL2")
  }
  end <- substr(tmp[2], start =0, stop = nchar(tmp[2])-1)
  str2 <- paste0(start, "_", end, "-1")
  str2
})
colnames(spliced_mat) <- cellID2
colnames(unspliced_mat) <- cellID2
rownames(spliced_mat) <- geneName
rownames(unspliced_mat) <- geneName
rownames(umap_mat) <- cellID2
names(clusters) <- cellID2

# format the MAESTRO matrices
colnames(maestro_day0) <- paste0("day0_", colnames(maestro_day0))
colnames(maestro_day10) <- paste0("day10_COCL2_", colnames(maestro_day10))
colnames(maestro_week5) <- paste0("week5_COCL2_", colnames(maestro_week5))
maestro <- cbind(maestro_day0, maestro_day10, maestro_week5)

# make sure all the cells are the same
final_cell_names <- intersect(colnames(maestro), intersect(colnames(all_data), colnames(spliced_mat)))
keep_vec <- rep(0, ncol(all_data_COCL2))
keep_vec[which(colnames(all_data_COCL2) %in% final_cell_names)] <- 1
all_data_COCL2$keep <- keep_vec
all_data_COCL2 <- subset(all_data_COCL2, keep == 1)

spliced_mat <- spliced_mat[,colnames(all_data_COCL2)]
unspliced_mat <- unspliced_mat[,colnames(all_data_COCL2)]
maestro <- maestro[,colnames(all_data_COCL2)]
umap_mat <- umap_mat[colnames(all_data_COCL2),]
clusters <- clusters[colnames(all_data_COCL2)]

all_data_COCL2[["spliced"]] <- Seurat::CreateAssayObject(data = spliced_mat)
all_data_COCL2[["spliced"]]@var.features <- intersect(rownames(spliced_mat), all_data_COCL2[["RNA"]]@var.features)
all_data_COCL2[["unspliced"]] <- Seurat::CreateAssayObject(data = unspliced_mat)
all_data_COCL2[["unspliced"]]@var.features <- intersect(rownames(unspliced_mat), all_data_COCL2[["RNA"]]@var.features)
all_data_COCL2[["maestro"]] <- Seurat::CreateAssayObject(data = maestro)
all_data_COCL2[["maestro"]]@var.features <- intersect(rownames(maestro), all_data_COCL2[["RNA"]]@var.features)
all_data_COCL2[["scVelo.umap"]] <- Seurat::CreateDimReducObject(umap_mat, key = "scVeloUMAP")
all_data_COCL2$scVelo_clusters <- clusters

save(all_data_COCL2, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_COCL2.RData")





