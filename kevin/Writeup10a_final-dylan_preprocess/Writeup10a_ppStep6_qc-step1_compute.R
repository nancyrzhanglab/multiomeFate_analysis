rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep5_barcode-assignment.RData"))

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

print("Set the variable genes in RNA")
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

print("Basic processing of RNA")
all_data <- Seurat::NormalizeData(all_data)
all_data <- Seurat::CellCycleScoring(all_data, 
                                     g2m.features = cc.genes$g2m.genes, 
                                     s.features = cc.genes$s.genes)
all_data[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = all_data, pattern = "^MT-")
all_data[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = all_data, pattern = "^RPS")

print("Basic processing of ATAC")
file_prefix <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/atac_fragments.tsv.gz"
file_folders <- c("2022_05_19_arc_time0", "2022_05_19_arc_time10_CIS", 
                  "2022_05_19_arc_time10_COCL2", "2022_05_19_arc_time10_DABTRAM",
                  "2022_05_19_arc_week5_CIS", "2022_05_19_arc_week5_COCL2",
                  "2022_05_19_arc_week5_DABTRAM")

Seurat::DefaultAssay(all_data) <- "ATAC"

# updating the path
for(i in 1:length(all_data[["ATAC"]]@fragments)){
  idx <- which(sapply(file_folders, function(x){length(grep(x, all_data[["ATAC"]]@fragments[[i]]@path)) > 0}))
  all_data[["ATAC"]]@fragments[[i]]@path <- paste0(file_prefix, file_folders[idx], file_suffix)
}

all_data <- Signac::NucleosomeSignal(object = all_data)
all_data <- Signac::TSSEnrichment(all_data, fast = FALSE)
all_data$high.tss <- ifelse(all_data$TSS.enrichment > 2, 'High', 'Low')
all_data$blacklist_fraction <- Signac::FractionCountsInRegion(
  object = all_data, 
  assay = "ATAC",
  regions = blacklist_hg19
)

print("Finished preprocessing data")
save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep6_qc.RData"))
