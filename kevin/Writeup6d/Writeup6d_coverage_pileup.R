rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("coverage_extractor_singlecell.R")
source("coverage_extractor_singlecell-plotter.R")

source("../Writeup6b/gene_list.R")
source("gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))
gene_vec <- gene_vec[-73] # TEMPORARY 
