rm(list=ls())
library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0_chromvar_lightweight_noATAC.RData")

Seurat::DefaultAssay(all_data) <- "chromvar"
all_data[["chromvar"]]@data[1:5,1:5]
dim(all_data[["chromvar"]]@data)