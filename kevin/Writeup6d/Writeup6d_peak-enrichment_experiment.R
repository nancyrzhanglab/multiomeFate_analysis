rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

quantile(all_data$TSS.percentile)
zz <- all_data[["ATAC"]]
zz2 <- zz@positionEnrichment
class(zz2)
names(zz2)
dim(zz2[["TSS"]])
class(zz2[["TSS"]])
zz2[["TSS"]][1:5,1:5]
zz2[["TSS"]][(nrow(zz2[["TSS"]])-5):nrow(zz2[["TSS"]]),1:5]
head(colnames(zz2[["TSS"]]))
tail(colnames(zz2[["TSS"]]))
Signac::Fragments(all_data[["ATAC"]])
