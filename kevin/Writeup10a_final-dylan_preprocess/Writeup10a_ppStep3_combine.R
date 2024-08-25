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

print("Adjusting the RNA dataset")
merged_counts <- do.call(cbind, lapply(1:7, function(i){
  SeuratObject::LayerData(all_data_rna, 
                          assay = "RNA",
                          layer = paste0("counts.", i))
}))
SeuratObject::LayerData(all_data_rna, 
                        assay = "RNA",
                        layer = "counts") <- merged_counts

for(i in 1:7){
  SeuratObject::LayerData(all_data_rna, 
                          assay = "RNA",
                          layer = paste0("counts.", i)) <- NULL
  SeuratObject::LayerData(all_data_rna, 
                          assay = "RNA",
                          layer = paste0("data.", i)) <- NULL
}

print("ATAC")
print(all_data_atac)
print("RNA")
print(all_data_rna)

atac_cells <- Seurat::Cells(all_data_atac)
rna_cells <- Seurat::Cells(all_data_rna)
common_cells <- intersect(atac_cells, rna_cells)

all_data_atac <- subset(all_data_atac, cells = common_cells)
all_data_rna <- subset(all_data_rna, cells = common_cells)
gc(TRUE)

print("ATAC")
print(all_data_atac)
print("RNA")
print(all_data_rna)

atac_cells <- Seurat::Cells(all_data_atac)
rna_cells <- Seurat::Cells(all_data_rna)
print("ATAC")
print(head(atac_cells))
print("RNA")
print(head(rna_cells))
all(atac_cells == rna_cells)

rna_counts <- SeuratObject::LayerData(all_data_rna, 
                                      assay = "RNA",
                                      layer = "counts")
rna_counts <- rna_counts[,Seurat::Cells(all_data_atac)]
all_data_atac[["RNA"]] <- Seurat::CreateAssayObject(rna_counts,
                                                    check.matrix = TRUE)
all_data <- all_data_atac
rm(list = c("all_data_atac", "all_data_rna")); gc(TRUE)

print(table(all_data$dataset))
print(head(all_data@meta.data))

print("Saving")
save(all_data, date_of_run, session_info,
     file = paste0(out_folder, "Writeup10a_ppStep3_combine.RData"))

########

print("RNA")
zz <- SeuratObject::LayerData(all_data, 
                              assay = "RNA",
                              layer = "counts")
print(zz[1:5,1:5])
print(quantile(zz@x))

print("ATAC")
zz <- SeuratObject::LayerData(all_data, 
                              assay = "ATAC",
                              layer = "counts")
print(zz[1:5,1:5])
print(quantile(zz@x))


