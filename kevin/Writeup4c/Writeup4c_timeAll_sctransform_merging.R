rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

create_seurat_object <- function(file_prefix, file_folder, file_suffix){
  tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  fragpath <- paste0(file_prefix, file_folder, "/outs/atac_fragments.tsv.gz")
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = tmp[["Gene Expression"]],
    assay = "RNA"
  )
  
  seurat_obj <- Seurat::SCTransform(seurat_obj, 
                                    method = "glmGamPoi", 
                                    variable.features.n = 500,
                                    verbose = T)
  
  seurat_obj
}

########

file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_04_arc_time0", "2022_04_arc_time10_CIS", 
                  "2022_04_arc_time10_COCL2", "2022_04_arc_time10_DABTRAM",
                  "2022_04_arc_week5_CIS", "2022_04_arc_week5_COCL2",
                  "2022_04_arc_week5_DABTRAM")
name_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM",
              "week5_CIS", "week5_COCL2", "week5_DABTRAM")

seurat_list <- lapply(1:length(file_folders), function(i){
  print(i)
  file_folder <- file_folders[i]
  seurat_obj <- create_seurat_object(file_prefix = file_prefix,
                                     file_folder = file_folder,
                                     file_suffix = file_suffix)
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})

# determine the important genes
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
gene_list <- lapply(seurat_list, function(seurat_obj){
  seurat_obj[["SCT"]]@var.features
})
gene_vec <- sort(unique(c(jackpot_genes, unlist(gene_list))))
all_genes <- unique(unlist(lapply(seurat_list, function(seurat_obj){
  rownames(seurat_obj[["RNA"]])
})))
gene_vec <- intersect(gene_vec, all_genes)

# run SCTransform separately on each dataset
seurat_list2 <- lapply(seurat_list, function(seurat_obj){
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj[["SCT"]] <- NULL
  
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = seurat_obj, pattern = "^RPS")
  
  seurat_obj <- Seurat::SCTransform(seurat_obj, 
                                    method = "glmGamPoi", 
                                    residual.features = gene_vec,
                                    verbose = T)
  seurat_obj
})
sapply(seurat_list2, function(seurat_obj){print(dim(seurat_obj[["SCT"]]@scale.data))})
sapply(seurat_list2, function(seurat_obj){print(length(seurat_obj[["SCT"]]@var.features))})

# genes that have almost all 0's are omitted, so we need to add those genes back
# into scale.data prior to merging
seurat_list3 <- lapply(seurat_list2, function(seurat_obj){
  mat <- matrix(0, 
                ncol = ncol(seurat_obj[["SCT"]]@scale.data), 
                nrow = length(gene_vec))
  colnames(mat) <- colnames(seurat_obj[["SCT"]]@scale.data)
  rownames(mat) <- gene_vec
  mat[rownames(seurat_obj[["SCT"]]@scale.data),] <- seurat_obj[["SCT"]]@scale.data
  seurat_obj[["SCT"]]@scale.data <- mat
  seurat_obj
})

all_data <- merge(seurat_list3[[1]], y = c(seurat_list3[[2]], seurat_list3[[3]], 
                                           seurat_list3[[4]], seurat_list3[[5]], 
                                           seurat_list3[[6]], seurat_list3[[7]]), 
                  add.cell.ids = name_vec, 
                  project = "All_Data", merge.data = T)
all_data[["SCT"]]@var.features <- gene_vec
dim(all_data[["SCT"]]@scale.data)
dim(all_data[["SCT"]]@data)

Seurat::DefaultAssay(all_data) <- "SCT"
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE) 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 1:50)

all_data <- Seurat::CellCycleScoring(all_data, 
                                     g2m.features = cc.genes$g2m.genes, 
                                     s.features = cc.genes$s.genes)

Seurat::DefaultAssay(all_data) <- "SCT"
plot1 <-Seurat::DimPlot(all_data, reduction = "umap",
                        group.by = "dataset", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("RNA (SCTransform),\n", length(all_data[["SCT"]]@var.features), " genes, using 50 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4c/Writeup4c_rna_umap_sctransform.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(all_data, 
                             reduction = "umap",
                             features = sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E", 
                                               "CD44", "LOXL2", "ID3")),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4c/Writeup4c_rna_umap_sctransform_jackpot1.png"),
                plot1, device = "png", width = 12, height = 12, units = "in")

save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4c/Writeup4c_timeAll_sctransform_merging.RData")


