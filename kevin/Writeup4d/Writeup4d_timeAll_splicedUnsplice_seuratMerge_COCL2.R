rm(list=ls())

library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_peakmerging.RData")

load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_fasttopics_COCL2.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_timeAll_saver_onlyRNA_COCL2-subset.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_spliced_mat_COCL2_full.RData")
load("../../../../out/kevin/Writeup4d/Writeup4d_unspliced_mat_COCL2_full.RData")
maestro_day0 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time0_202204/MAESTRO_arc0_rpmatrix.RDS")
maestro_day10 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/time10_cocl2_202204/MAESTRO_arc10_cocl2_rpmatrix.RDS")
maestro_week5 <- readRDS("~/project/Multiome_fate/BarcodeOutputs/2022_02/Sijia_outputs/maestro/week5_cocl2_202204/MAESTRO_week5_cocl2_rpmatrix.RDS")

load("../../../../out/kevin/Writeup4d/Writeup4d_scvelo_umap_COCL2.RData")
cellIDb <- cellID
geneNameb <- geneName

load("../../../../out/kevin/Writeup4d/Writeup4d_info_COCL2_full.RData")
stopifnot(all(cellIDb == cellID))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

keep_vec <- rep(0, ncol(all_data))
keep_vec[which(colnames(all_data) %in% colnames(all_data_subset))] <- 1
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == 1)
all_data_subset[["ATAC"]] <- all_data[["ATAC"]]

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
final_cell_names <- intersect(colnames(maestro), intersect(colnames(all_data_subset), colnames(spliced_mat)))
keep_vec <- rep(0, ncol(all_data_subset))
keep_vec[which(colnames(all_data_subset) %in% final_cell_names)] <- 1
all_data_subset$keep <- keep_vec
all_data_subset <- subset(all_data_subset, keep == 1)

spliced_mat <- spliced_mat[,colnames(all_data_subset)]
unspliced_mat <- unspliced_mat[,colnames(all_data_subset)]
maestro <- maestro[,colnames(all_data_subset)]
umap_mat <- umap_mat[colnames(all_data_subset),]
clusters <- clusters[colnames(all_data_subset)]

all_data_subset[["spliced"]] <- Seurat::CreateAssayObject(counts = spliced_mat)
all_data_subset[["spliced"]]@var.features <- intersect(rownames(spliced_mat), all_data_subset[["RNA"]]@var.features)
all_data_subset[["unspliced"]] <- Seurat::CreateAssayObject(counts = unspliced_mat)
all_data_subset[["unspliced"]]@var.features <- intersect(rownames(unspliced_mat), all_data_subset[["RNA"]]@var.features)
all_data_subset[["maestro"]] <- Seurat::CreateAssayObject(counts = maestro)
all_data_subset[["maestro"]]@var.features <- intersect(rownames(maestro), all_data_subset[["RNA"]]@var.features)
all_data_subset[["scVelo.umap"]] <- Seurat::CreateDimReducObject(umap_mat, key = "scVeloUMAP")
all_data_subset$scVelo_clusters <- clusters

topic_res$L <- topic_res$L[colnames(all_data_subset),]
all_data_subset[["fasttopic"]] <- Seurat::CreateDimReducObject(embeddings = topic_res$L, 
                                                               loadings =  topic_res$F,
                                                               assay = "RNA",
                                                               key = "fastTopic_")

Seurat::DefaultAssay(all_data_subset) <- "ATAC"
gene_activities <- Signac::GeneActivity(all_data_subset)
all_data_subset[["geneActivity"]] <- CreateAssayObject(counts = gene_activities)
all_data_subset[["geneActivity"]]@var.features <- intersect(rownames(gene_activities), all_data_subset[["RNA"]]@var.features)
all_data_subset <- Seurat::NormalizeData(
  object = all_data_subset,
  assay = "geneActivity",
  normalization.method = "LogNormalize",
  scale.factor = median(all_data_subset$nCount_geneActivity)
)

##########3

print("Seurat preprocess of ATAC")
Seurat::DefaultAssay(all_data_subset) <- "ATAC"
set.seed(10)
all_data_subset <- Signac::RunTFIDF(all_data_subset)
all_data_subset <- Signac::FindTopFeatures(all_data_subset, min.cutoff = 'q0')
all_data_subset <- Signac::RunSVD(all_data_subset)
set.seed(10)
all_data_subset <- Seurat::RunUMAP(object = all_data_subset, 
                                   reduction = 'lsi', dims = 2:50,
                                   reduction.name = 'atac.umap')

# https://satijalab.org/seurat/archive/v3.0/atacseq_integration_vignette.html
print("Seurat preprocess of gene activity")
Seurat::DefaultAssay(all_data_subset) <- "geneActivity"
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)
all_data_subset <- Seurat::RunPCA(all_data_subset, verbose = FALSE,
                                  reduction.name = "activityPCA") 
set.seed(10)
all_data_subset <- Seurat::RunUMAP(all_data_subset, dims = 2:50,
                                   reduction = "activityPCA", 
                                   reduction.name = "activity.umap")

print("Seurat preprocess of Maestro")
Seurat::DefaultAssay(all_data_subset) <- "maestro"
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)
all_data_subset <- Seurat::RunPCA(all_data_subset, verbose = FALSE,
                                  reduction.name = "maestroPCA") 
set.seed(10)
all_data_subset <- Seurat::RunUMAP(all_data_subset, dims = 2:50,
                                   reduction = "maestroPCA", 
                                   reduction.name = "maestro.umap")


### update the RNA things too
Seurat::DefaultAssay(all_data_subset) <- "RNA"
all_data_subset <- Seurat::NormalizeData(all_data_subset)
all_data_subset <- Seurat::ScaleData(all_data_subset)
all_data_subset <- Seurat::RunPCA(all_data_subset, verbose = FALSE) 
set.seed(10)
all_data_subset <- Seurat::RunUMAP(all_data_subset, dims = 1:50)


Seurat::DefaultAssay(all_data_subset) <- "Saver"
all_data_subset <- Seurat::ScaleData(all_data_subset)
all_data_subset <- Seurat::RunPCA(all_data_subset, verbose = FALSE,
                                  reduction.name = "saverpca") 
set.seed(10)
all_data_subset <- Seurat::RunUMAP(all_data_subset, dims = 1:50,
                                   reduction = "saverpca",
                                   reduction.name = "saverumap")

##### normalize spliced/unspliced
Seurat::DefaultAssay(all_data_subset) <- "Saver"
all_data_subset <- Seurat::FindNeighbors(all_data_subset, 
                                         reduction = "saverpca",
                                         dims = 1:50)

depth <- Matrix::colSums(all_data_subset[["RNA"]]@counts)
depth <- pmax(depth, 1)
median_depth <- stats::median(depth)
genes <- all_data_subset[["spliced"]]@var.features
mat <- all_data_subset[["spliced"]]@counts
mat <- tiltedCCA:::.mult_mat_vec(mat, median_depth/depth)
mat <- mat %*% all_data_subset@graphs[["Saver_snn"]] 
all_data_subset[["spliced"]]@data <- mat

genes <- all_data_subset[["unspliced"]]@var.features
mat <- all_data_subset[["unspliced"]]@counts
mat <- tiltedCCA:::.mult_mat_vec(mat, median_depth/depth)
mat <- mat %*% all_data_subset@graphs[["Saver_snn"]] 
all_data_subset[["unspliced"]]@data <- mat

save(all_data_subset, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4d/Writeup4d_timeAll_splicedUnspliced_seuratMerge_COCL2.RData")

#########

plot1 <-Seurat::DimPlot(all_data_subset, reduction = "atac.umap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2: ATAC"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_atac_umap-cocl2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-Seurat::DimPlot(all_data_subset, reduction = "maestro.umap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2: ATAC (Maestro)"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_atac-maestro_umap-cocl2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-Seurat::DimPlot(all_data_subset, reduction = "activity.umap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2: ATAC (Gene activity)"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_atac-activity_umap-cocl2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-Seurat::DimPlot(all_data_subset, reduction = "umap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2: RNA (Log-normalized)"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_rna_umap-cocl2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-Seurat::DimPlot(all_data_subset, reduction = "saverumap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("COCL2: RNA (SAVER)"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4d/Writeup4d_rna-saver_umap-cocl2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

