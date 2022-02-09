rm(list=ls())
load("../../../../data/Sydney_stressors/2021-10-28/2021_10_28_Cleaned_Data/final_all_data.RData")
load("../../../../data/Sydney_stressors/2021-10-28/2021_10_28_Cleaned_Data/final_all_data2.RData")

###########

lineage_names <- rownames(all_data2[["lineage"]]@counts)
zz <- all_data[["lineage"]]@counts
zz <- zz[rownames(zz) %in% lineage_names,]

cell_names <- rownames(all_data@meta.data)
naive_cells <- cell_names[which(all_data$OG_condition=="naive")]
naive_idx <- which(colnames(zz) %in% naive_cells)
which(sparseMatrixStats::rowSums2(zz[,naive_idx]) == 0)