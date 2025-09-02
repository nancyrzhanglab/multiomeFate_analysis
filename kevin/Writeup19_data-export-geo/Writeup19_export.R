rm(list=ls())
library(Seurat)
library(Signac)
library(multiomeFate)
library(Matrix)

out_folder <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/out/kevin/Writeup19/"

# all_data <- multiomeFate:::data_loader(which_files = c("atac", "lineage", "rna", "fasttopics", "peakvi", "rna_dimred", "wnn", "fatepotential"))


all_data <- multiomeFate:::data_loader(which_files = c("atac"))
# https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/externalFormats.html
count_mat <- SeuratObject::LayerData(all_data,
                                     layer = "counts",
                                     assay = "ATAC")
write.table(colnames(count_mat), 
            file = paste0(out_folder, "cells_atac.tsv"),
            quote = FALSE, 
            sep = '\t', 
            row.names = FALSE,
            col.names = FALSE)
write.table(rownames(count_mat), 
            file = paste0(out_folder, "features_atac.tsv"),
            quote = FALSE, 
            sep = '\t', 
            row.names = FALSE,
            col.names = FALSE)
Matrix::writeMM(count_mat, 
                file = paste0(out_folder, "counts_atac.mtx"))

df <- as.data.frame(all_data[["ATAC"]]@ranges)
rm(ls = "all_data"); gc(TRUE)

all_data <- multiomeFate:::data_loader(which_files = c("lineage", "rna"))
count_mat <- SeuratObject::LayerData(all_data,
                                     layer = "counts",
                                     assay = "Lineage")
write.table(colnames(count_mat), 
            file = paste0(out_folder, "cells_lineage.tsv"),
            quote = FALSE, 
            sep = '\t', 
            row.names = FALSE,
            col.names = FALSE)
write.table(rownames(count_mat), 
            file = paste0(out_folder, "features_lineage.tsv"),
            quote = FALSE, 
            sep = '\t', 
            row.names = FALSE,
            col.names = FALSE)
Matrix::writeMM(count_mat, 
                file = paste0(out_folder, "counts_lineage.mtx"))

count_mat <- SeuratObject::LayerData(all_data,
                                     layer = "counts",
                                     assay = "RNA")
write.table(colnames(count_mat), 
            file = paste0(out_folder, "cells_rna.tsv"),
            quote = FALSE, 
            sep = '\t', 
            row.names = FALSE,
            col.names = FALSE)
write.table(rownames(count_mat), 
            file = paste0(out_folder, "features_rna.tsv"),
            quote = FALSE, 
            sep = '\t', 
            row.names = FALSE,
            col.names = FALSE)
Matrix::writeMM(count_mat, 
                file = paste0(out_folder, "counts_rna.mtx"))
rm(ls = "all_data"); gc(TRUE)

all_data <- multiomeFate:::data_loader(which_files = c("saver", "fasttopics", "peakvi", "rna_dimred", "wnn", "fatepotential"))

fasttopics_cis <- all_data[["fasttopic.CIS"]]@cell.embeddings
fasttopics_cocl2 <- all_data[["fasttopic.COCL2"]]@cell.embeddings
fasttopics_dabtram <- all_data[["fasttopic.DABTRAM"]]@cell.embeddings
fasttopics <- cbind(fasttopics_cis, fasttopics_cocl2, fasttopics_dabtram)
write.csv(fasttopics,
          file = paste0(out_folder, "dimred_fasttopics.csv"))

peakvi_cis <- all_data[["peakVI.CIS"]]@cell.embeddings
peakvi_cocl2 <- all_data[["peakVI.COCL2"]]@cell.embeddings
peakvi_dabtram <- all_data[["peakVI.DABTRAM"]]@cell.embeddings

peakvi <- matrix(NA, 
                 nrow = length(Seurat::Cells(all_data)),
                 ncol = ncol(peakvi_cis) + ncol(peakvi_cocl2) + ncol(peakvi_dabtram),
                 dimnames = list(Seurat::Cells(all_data), 
                                 c(colnames(peakvi_cis), colnames(peakvi_cocl2), colnames(peakvi_dabtram))))
peakvi[rownames(peakvi_cis), colnames(peakvi_cis)] <- peakvi_cis
peakvi[rownames(peakvi_cocl2), colnames(peakvi_cocl2)] <- peakvi_cocl2
peakvi[rownames(peakvi_dabtram), colnames(peakvi_dabtram)] <- peakvi_dabtram
write.csv(peakvi,
          file = paste0(out_folder, "dimred_peakvi.csv"))

umap <- all_data[["Saver.umap"]]@cell.embeddings
write.csv(umap,
          file = paste0(out_folder, "umap_rna.csv"))
umap <- all_data[["pVI.All.umap"]]@cell.embeddings
write.csv(umap,
          file = paste0(out_folder, "umap_atac.csv"))
umap <- all_data[["wnn.umap"]]@cell.embeddings
write.csv(umap,
          file = paste0(out_folder, "umap_wnn.csv"))


metadata <- all_data@meta.data
write.csv(metadata,
          file = paste0(out_folder, "metadata.csv"))



