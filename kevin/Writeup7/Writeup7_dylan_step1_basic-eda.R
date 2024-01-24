rm(list=ls())
library(Seurat)

load("~/nzhanglab/data/DylanSchaff/all_data_lineages.RData")
set.seed(10)

print("Remove certain cells")
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[intersect(which(!is.na(all_data$assigned_lineage)),
                   which(all_data$assigned_posterior >= 0.5))] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

print("Seurat preprocess of RNA")
Seurat::DefaultAssay(all_data) <- "RNA"
treatment_vec <- sort(unique(all_data$OG_condition))
  
var_list <- lapply(treatment_vec, function(treatment){
  print(treatment)
  tmp <- all_data
  keep_vec <- rep(FALSE, ncol(tmp))
  keep_vec[tmp$OG_condition %in% treatment] <- TRUE
  tmp$keep <- keep_vec
  tmp <- subset(tmp, keep == 1)
  
  set.seed(10)
  tmp <- Seurat::FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 500)
  vec1 <- Seurat::VariableFeatures(tmp) 
  
  if(treatment != "naive") {
    tmp <- all_data
    keep_vec <- rep(FALSE, ncol(tmp))
    keep_vec[tmp$OG_condition %in% c(treatment, "naive")] <- TRUE
    tmp$keep <- keep_vec
    tmp <- subset(tmp, keep == 1)
    
    set.seed(10)
    tmp <- Seurat::FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 500)
    vec2 <- Seurat::VariableFeatures(tmp) 
  } else {
    vec2 <- numeric(0)
  }
  
  unique(c(vec1, vec2))
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
Seurat::VariableFeatures(all_data) <- all_genes

save(all_data, 
     file = "~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step1_basic-eda.RData")

all_data <- Seurat::NormalizeData(all_data)
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE) 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 1:30)

print("Finished preprocessing data")
date_of_run <- Sys.time()
session_info <- devtools::session_info()
save(all_data, date_of_run, session_info,
     file = "~/project/Multiome_fate/out/kevin/Writeup7/Writeup7_dylan_step1_basic-eda.RData")

print("Done! :)")