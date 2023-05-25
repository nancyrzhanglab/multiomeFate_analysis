rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

all_data <- subset(all_data, dataset %in% c("day0", "day10_DABTRAM", "week5_DABTRAM"))
all_data <- subset(all_data, assigned_posterior >= 0.5)
all_data2 <- subset(all_data, dataset == "day10_DABTRAM")

###################

# do a pseudotime analysis just on the Day10s 
