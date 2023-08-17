rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6b/Writeup6b_chromVar.RData")

all_data

print("Simplifying dataset")
Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["geneActivity"]] <- NULL

all_data[["pca"]] <- NULL
all_data[["umap"]] <- NULL
all_data[["lsi"]] <- NULL
all_data[["atac.umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL

print("Removing cells without a barcode")
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[intersect(which(!is.na(all_data$assigned_lineage)),
                   which(all_data$assigned_posterior >= 0.5))] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

print("Keeping only day0 cells")
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[all_data$dataset == "day0"] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

print("Saving")
save(date_of_run, session_info, 
     all_data, 
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract.RData")