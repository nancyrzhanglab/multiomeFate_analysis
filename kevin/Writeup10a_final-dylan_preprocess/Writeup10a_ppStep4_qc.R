
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

all_data <- Seurat::NormalizeData(all_data)
all_data <- Seurat::CellCycleScoring(all_data, 
                                     g2m.features = cc.genes$g2m.genes, 
                                     s.features = cc.genes$s.genes)

print("Finished preprocessing data")
save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4e/Writeup4e_timeAll_peakmerging.RData")

