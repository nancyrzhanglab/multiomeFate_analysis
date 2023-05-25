rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6h/Writeup6h_DABTRAM-day10_pseudotime.RData")
all_data_dabtram <- all_data
load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

# where are the day0 cells that do not expand, and how many unique lineages do they belong to?

# do the chromvar analysis 

Seurat::DefaultAssay(all_data) <- "chromvar"
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = "winning_day0",
  ident.2 = "losing_day0",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  verbose = F
)