rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# remove all the ATAC things
Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["lsi"]] <- NULL
all_data[["geneActivity"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

treatment_vec <- c("CIS", "COCL2", "DABTRAM")
for(treatment in treatment_vec){
  print(treatment)
  peakvi_mat <- read.csv(paste0("../../../../out/kevin/Writeup6o/Writeup6o_all-data-atac_", treatment, "_peakVI.csv"),
                         row.names = 1)
  svd_res <- svd(peakvi_mat)
  tmp <- sweep(svd_res$u, MARGIN = 2, STATS = svd_res$d, FUN = "*")
  rownames(tmp) <- rownames(peakvi_mat)
  peakvi_mat <- tmp
  peakvi_mat <- scale(peakvi_mat)
  colnames(peakvi_mat) <- paste0("peakVI", treatment, "_", 1:ncol(peakvi_mat))
  
  # create an empty dimred
  dimred_mat <- matrix(NA, nrow = ncol(all_data), ncol = ncol(peakvi_mat))
  rownames(dimred_mat) <- colnames(all_data)
  colnames(dimred_mat) <- paste0("peakVI", treatment, "_", 1:ncol(dimred_mat))
  
  # put in the embedding
  dimred_mat[rownames(peakvi_mat),] <- peakvi_mat
  
  all_data[[paste0("peakVI", treatment)]] <- Seurat::CreateDimReducObject(dimred_mat,
                                                                          assay = "ATAC")
}

print("Saving")
save(date_of_run, session_info,
     all_data,
     file = "../../../../out/kevin/Writeup6o/Writeup6o_all-data_peakvi.RData")

print("Done! :)")
