rm(list=ls())
library(Seurat)

load("../../../../out/kevin/Writeup6p/Writeup6p_all-data_lightweight_noATAC.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

treatment <- "COCL2"
day_early <- "day10"
day_later <- "week5"
day_early_full <- paste0(day_early, "_", treatment)
day_later_full <- paste0(day_later, "_", treatment)

# keep only the relevant cells
keep_vec <- rep(FALSE, ncol(all_data))
idx <- which(all_data$dataset %in% c("day0", paste0("day10_", treatment), paste0("week5_", treatment)))
keep_vec[idx] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

set.seed(10)
all_data <- Seurat::RunUMAP(all_data, 
                            reduction = paste0("fasttopic_", treatment),
                            dims = 1:30)

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["Saver"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["peakVI_CIS"]] <- NULL
all_data[["peakVI_DABTRAM"]] <- NULL
all_data[["pca"]] <- NULL
all_data[["atac.umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["saverumap"]] <- NULL

all_data[["umap"]]@assay.used <- "Lineage"
all_data[["fasttopic_COCL2"]]@assay.used <- "Lineage"
all_data[["peakVI_COCL2"]]@assay.used <- "Lineage"

Seurat::DefaultAssay(all_data) <- "Lineage"

note <- "Lightweight data only for simulation purposes. Kept the RNA modality"
save(all_data, date_of_run, session_info, note,
     file = "~/project/Multiome_fate/out/kevin/Writeup8b/Writeup8b_simulation_day10-COCL2_withRNA.RData")




