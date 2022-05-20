rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_peakmerging.RData")

maestro_day0 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time0_202204/MAESTRO_arc0_rpmatrix.RDS")
maestro_day10 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time10_cis_202204/MAESTRO_arc10_cis_rpmatrix.RDS")
maestro_week5 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/week5_cis_202204/MAESTRO_week5_cis_rpmatrix.RDS")

# format the MAESTRO matrices
colnames(maestro_day0) <- paste0("day0_", colnames(maestro_day0))
colnames(maestro_day10) <- paste0("day10_CIS_", colnames(maestro_day10))
colnames(maestro_week5) <- paste0("week5_CIS_", colnames(maestro_week5))
maestro <- cbind(maestro_day0, maestro_day10, maestro_week5)

# make sure all the cells are the same
final_cell_names <- intersect(colnames(maestro), colnames(all_data))
keep_vec <- rep(0, ncol(all_data))
keep_vec[which(colnames(all_data) %in% final_cell_names)] <- 1
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == 1)

maestro <- maestro[,colnames(all_data)]
all_data[["maestro"]] <- Seurat::CreateAssayObject(counts = maestro)
all_data[["maestro"]]@var.features <- intersect(rownames(maestro), all_data[["RNA"]]@var.features)

# see https://satijalab.org/signac/articles/pbmc_vignette.html
Seurat::DefaultAssay(all_data) <- "ATAC"
gene_activities <- Signac::GeneActivity(all_data)
all_data[["geneActivity"]] <- CreateAssayObject(counts = gene_activities)
all_data <- Seurat::NormalizeData(
  object = all_data,
  assay = "geneActivity",
  normalization.method = "LogNormalize",
  scale.factor = median(all_data$nCount_geneActivity)
)

########################

print("Seurat preprocess of ATAC")
Seurat::DefaultAssay(all_data) <- "ATAC"
set.seed(10)
all_data <- Signac::RunTFIDF(all_data)
all_data <- Signac::FindTopFeatures(all_data, min.cutoff = 'q0')
all_data <- Signac::RunSVD(all_data)
set.seed(10)
all_data <- Seurat::RunUMAP(object = all_data, 
                                   reduction = 'lsi', dims = 2:50,
                                   reduction.name = 'atac.umap')



