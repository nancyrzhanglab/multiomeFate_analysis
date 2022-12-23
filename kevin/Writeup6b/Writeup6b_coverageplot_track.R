rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

object <- all_data
region <- "ANXA1"
which_ident <- sort(unique(all_data$dataset))
assay = "ATAC"
cells = NULL
extend.downstream = 1000
extend.upstream = 1000
group.by = NULL
idents = NULL
multipler = 1e6
scale.factor = NULL
sep = c("-", "-")
window = 100

Seurat::Idents(object) <- "dataset"