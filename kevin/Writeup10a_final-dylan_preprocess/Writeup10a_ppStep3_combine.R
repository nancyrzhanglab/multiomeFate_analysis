rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

print("Combine both dataset")
out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep1_peakmerging.RData"))
load(paste0(out_folder, "Writeup10a_ppStep2_rna-merge.RData"))

print(all_data_atac)
print(all_data_rna)

atac_cells <- Seurat::Cells(all_data_atac)
rna_cells <- Seurat::Cells(all_data_rna)
common_cells <- intersect(atac_cells, rna_cells)

all_data_atac <- subset(all_data_atac, cells = common_cells)
all_data_rna <- subset(all_data_rna, cells = common_cells)
gc(TRUE)

print(all_data_atac)
print(all_data_rna)

atac_cells <- Seurat::Cells(all_data_atac)
rna_cells <- Seurat::Cells(all_data_rna)
print(head(atac_cells))
print(head(rna_cells))
all(atac_cells == rna_cells)

all_data_atac[["RNA"]] <- all_data_rna
all_data <- all_data_atac
rm(list = c("all_data_atac", "all_data_rna")); gc(TRUE)

print(table(all_data$dataset))
print(head(all_data@meta.data))

print("Saving")
save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep3_combine.RData"))




