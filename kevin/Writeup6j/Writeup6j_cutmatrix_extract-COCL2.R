rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6f/gene_vec.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "COCL2"

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["geneActivity"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL
all_data[["activity.umap"]] <- NULL

keep_vec <- rep(NA, ncol(all_data))
idx1 <- which(all_data$dataset %in% c("day0", 
                                      paste0("day10_", treatment), 
                                      paste0("week5_", treatment)))
idx2 <- which(all_data$assigned_posterior >= 0.5)
idx3 <- which(!is.na(all_data$assigned_lineage))
keep_vec[intersect(intersect(idx1, idx2), idx3)] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

gene_vec_all <- sort(intersect(gene_vec_all,
                               rownames(all_data[["Saver"]]@scale.data)))

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

surviving_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
dying_lineages <- rownames(tab_mat)[intersect(
  which(tab_mat[,paste0("day10_", treatment)] <= 2),
  which(tab_mat[,paste0("week5_", treatment)] <= 2)
)]
winning_idx <- intersect(
  which(all_data$assigned_lineage %in% surviving_lineages),
  which(all_data$dataset == "day0")
)
dying_idx <- intersect(
  which(all_data$assigned_lineage %in% dying_lineages),
  which(all_data$dataset == "day0")
)
winning_cells <- colnames(all_data)[winning_idx]
dying_cells <- colnames(all_data)[dying_idx]

##############

cutmat_list <- vector("list", length(gene_vec_all))
names(cutmat_list) <- gene_vec_all

for(zz in 1:length(gene_vec_all)){
  if(zz %% 100 == 0) {
    save(cutmat_list, date_of_run, session_info,
         file = paste0("../../../../out/kevin/Writeup6i/day0_cutmatrix_extract-", treatment, "_tmp.RData"))
  }
  
  print(paste0(zz, " of ", length(gene_vec_all), " in ", treatment))
  gene <- gene_vec_all[zz]
  
  cutmat_winning <- multiomeFate:::extract_cutmatrix(
    object = all_data,
    gene = gene,
    cells = winning_cells
  )
  cutmat_dying <- multiomeFate:::extract_cutmatrix(
    object = all_data,
    gene = gene,
    cells = dying_cells
  )
  if(all(is.null(cutmat_winning)) || all(is.null(cutmat_dying))) next()
  cutmat_all <- rbind(cutmat_winning, cutmat_dying)
  
  peak_mat <- multiomeFate:::extract_peaks(
    object = all_data,
    gene = gene
  )
  if(nrow(peak_mat) == 0) next()
  
  peak_prior <-  multiomeFate:::compute_peak_prior(cutmat = cutmat_all,
                                                   peak_mat = peak_mat)
  if(all(is.na(peak_prior))) next()
  
  cutmat_res <- list(cutmat_winning = cutmat_winning,
                     cutmat_dying = cutmat_dying,
                     peak_mat = peak_mat,
                     peak_prior = peak_prior)
  cutmat_list[[gene]] <- cutmat_res
}

save(cutmat_list, date_of_run, session_info,
     file = paste0("../../../../out/kevin/Writeup6i/day0_cutmatrix_extract-", treatment, "_tmp.RData"))
