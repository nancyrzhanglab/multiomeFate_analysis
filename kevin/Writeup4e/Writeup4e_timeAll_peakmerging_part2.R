rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")

print("Seurat preprocess of RNA")
Seurat::DefaultAssay(all_data) <- "RNA"
subset_list <- list(day0 = "day0", 
                    day10_CIS = "day10_CIS",
                    day10_COCL2 = "day10_COCL2",
                    day10_DABTRAM = "day10_DABTRAM",
                    week5_CIS = "week5_CIS",
                    week5_COCL2 = "week5_COCL2",
                    week5_DABTRAM = "week5_DABTRAM",
                    CIS = c("day0", "day10_CIS", "week5_CIS"),
                    COCL2 = c("day0", "day10_COCL2", "week5_COCL2"),
                    DABTRAM = c("day0", "day10_DABTRAM", "week5_DABTRAM"),
                    all = c("day0", "day10_CIS", "week5_CIS",
                            "day10_COCL2", "week5_COCL2", "day10_DABTRAM", "week5_DABTRAM"))
var_list <- lapply(1:length(subset_list), function(i){
  print(i)
  tmp <- all_data
  tmp[["ATAC"]] <- NULL
  keep_vec <- rep(0, ncol(tmp))
  keep_vec[tmp$dataset %in% subset_list[[i]]] <- 1
  tmp$keep <- keep_vec
  tmp <- subset(tmp, keep == 1)
  
  set.seed(10)
  tmp <- Seurat::FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 500)
  tmp[["RNA"]]@var.features
})                   

jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3")
jackpot_genes <- intersect(jackpot_genes, rownames(all_data))
hk_genes <- read.csv("~/project/eSVD/data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
all_genes <- unique(c(unlist(var_list), jackpot_genes, hk_genes, cycling_genes))
all_genes <- intersect(all_genes, rownames(all_data))

all_data[["RNA"]]@var.features <- all_genes
save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")

all_data <- Seurat::NormalizeData(all_data)
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE) 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 1:50)
all_data <- Seurat::CellCycleScoring(all_data, 
                                     g2m.features = cc.genes$g2m.genes, 
                                     s.features = cc.genes$s.genes)

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

print("Seurat preprocess of gene activity")
Seurat::DefaultAssay(all_data) <- "ATAC"
gene_activities <- Signac::GeneActivity(all_data)
all_data[["geneActivity"]] <- CreateAssayObject(counts = gene_activities)
Seurat::DefaultAssay(all_data) <- "geneActivity"
all_data[["geneActivity"]]@var.features <- intersect(rownames(gene_activities), all_data[["RNA"]]@var.features)
all_data <- Seurat::NormalizeData(
  object = all_data,
  assay = "geneActivity",
  normalization.method = "LogNormalize",
  scale.factor = median(all_data$nCount_geneActivity)
)
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE,
                           reduction.name = "activityPCA") 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 2:50,
                            reduction = "activityPCA", 
                            reduction.name = "activity.umap")

print("Finished preprocessing data")
save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")

