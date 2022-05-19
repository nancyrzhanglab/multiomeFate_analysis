rm(list=ls())
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_peakmerging.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_scvelo_umap_CIS.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_spliced_mat_CIS.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_unspliced_mat_CIS.RData")
maestro_day0 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time0_202204/MAESTRO_arc0_rpmatrix.RDS")
maestro_day10 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time10_cis_202204/MAESTRO_arc10_cis_rpmatrix.RDS")
maestro_week5 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/week5_cis_202204/MAESTRO_week5_cis_rpmatrix.RDS")

library(Seurat)
library(Signac)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", "day10_CIS", "week5_CIS"))] <- 1
all_data$keep <- keep_vec
all_data_CIS <- subset(all_data, keep == 1)

#####################

# format the spliced and unspliced mat in
spliced_mat <- Matrix::Matrix(t(spliced_mat), sparse = T)
unspliced_mat <- Matrix::Matrix(t(unspliced_mat), sparse = T)
cellID2 <- sapply(cellID, function(str){
  tmp <- strsplit(str, split = ":")[[1]]
  tmp2 <- strsplit(tmp, split = "_")[[1]]
  start <- tmp2[3]
  if(start != "day0"){
    start <- paste0(start, "_CIS")
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
colnames(maestro_day10) <- paste0("day10_CIS_", colnames(maestro_day10))
colnames(maestro_week5) <- paste0("week5_CIS_", colnames(maestro_week5))
maestro <- cbind(maestro_day0, maestro_day10, maestro_week5)

# make sure all the cells are the same
final_cell_names <- intersect(colnames(maestro), intersect(colnames(all_data), colnames(spliced_mat)))
keep_vec <- rep(0, ncol(all_data_CIS))
keep_vec[which(colnames(all_data_CIS) %in% final_cell_names)] <- 1
all_data_CIS$keep <- keep_vec
all_data_CIS <- subset(all_data_CIS, keep == 1)

spliced_mat <- spliced_mat[,colnames(all_data_CIS)]
unspliced_mat <- unspliced_mat[,colnames(all_data_CIS)]
maestro <- maestro[,colnames(all_data_CIS)]
umap_mat <- umap_mat[colnames(all_data_CIS),]
clusters <- clusters[colnames(all_data_CIS)]

all_data_CIS[["spliced"]] <- Seurat::CreateAssayObject(data = spliced_mat)
all_data_CIS[["spliced"]]@var.features <- intersect(rownames(spliced_mat), all_data_CIS[["RNA"]]@var.features)
all_data_CIS[["unspliced"]] <- Seurat::CreateAssayObject(data = unspliced_mat)
all_data_CIS[["unspliced"]]@var.features <- intersect(rownames(unspliced_mat), all_data_CIS[["RNA"]]@var.features)
all_data_CIS[["maestro"]] <- Seurat::CreateAssayObject(data = maestro)
all_data_CIS[["maestro"]]@var.features <- intersect(rownames(maestro), all_data_CIS[["RNA"]]@var.features)
all_data_CIS[["scVelo.umap"]] <- Seurat::CreateDimReducObject(umap_mat, key = "scVeloUMAP")
all_data_CIS$scVelo_clusters <- clusters

save(all_data_CIS, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_CIS.RData")





