rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

tmp <- sapply(gene_vec, function(gene){
  tmp <- Signac::LookupGeneCoords(
    object = all_data,
    gene = gene,
    assay = "ATAC"
  )
  if(all(is.null(tmp))) return(F) else T
})
gene_vec <- gene_vec[which(tmp)]

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "COCL2"

# find the winning and losing cells
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
surviving_lineages <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
dying_lineages <- rownames(tab_mat)[which(apply(tab_mat,1,max)<=1)]
winning_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% surviving_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
dying_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% dying_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
length(winning_idx); length(dying_idx)
winning_cells <- colnames(all_data)[winning_idx]
dying_cells <- colnames(all_data)[dying_idx]

result_list <- vector("list", length = length(gene_vec))
names(result_list) <- gene_vec

###############################

for(i in 1:length(gene_vec)){
  gene <- gene_vec[i]
  print(gene)
  print("Computing cutmat")
  
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
  cutmat_all <- rbind(cutmat_winning, cutmat_dying)
  
  peak_mat <- multiomeFate:::extract_peaks(
    object = all_data,
    gene = gene
  )
  peak_locations <- multiomeFate:::compute_peak_locations(peak_mat)
  # peak_prior <-  multiomeFate:::compute_peak_prior(cutmat = cutmat_all,
  #                                                 peak_mat = peak_mat)
  peak_prior <-  compute_peak_prior(cutmat = cutmat_all,
                                    peak_mat = peak_mat,
                                    min_prior = 0)
  peak_width <- stats::median(apply(peak_mat, 1, diff))
  
  # remove peaks with less than a prior of 0.01
  if(any(peak_prior <= 0.01)){
    idx <- which(peak_prior <= 0.01)
    print(paste0("Removed " , length(idx), " peaks"))
    peak_locations <- peak_locations[-idx]
    peak_prior <- peak_prior[-idx]
    peak_prior <- peak_prior/sum(peak_prior)
  }
  
  print("Computing crossfit")
  set.seed(i)
  result_list[[i]] <- multiomeFate:::peak_testing(
    cutmat_dying = cutmat_dying, 
    cutmat_winning = cutmat_winning,
    peak_locations = peak_locations,
    max_iter = 25,
    peak_prior = peak_prior,
    peak_width = peak_width,
    verbose = 1
  )
  if(!all(is.na(result_list[[i]]))){
    print(paste0("P-value of ", round(min(1/result_list[[i]]$teststat,1),5)))
  } else {
    print("No computation")
  }
  
  save(result_list, date_of_run, session_info,
       file = "../../../../out/kevin/Writeup6f/COCL_entropy-crossfit_tmp.RData")
}