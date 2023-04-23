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

treatment <- "DABTRAM"

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
  bin_midpoints <- multiomeFate:::compute_bin_midpoints(peak_mat)
  bin_limits <- c(bin_midpoints[1] - abs(bin_midpoints[2]-bin_midpoints[1]),
                  bin_midpoints[7] + abs(bin_midpoints[2]-bin_midpoints[1]))
  peak_locations <- multiomeFate:::compute_peak_locations(peak_mat)
  peak_prior <- multiomeFate:::compute_peak_prior(mat = cutmat_all,
                                                  peak_mat = peak_mat)
  
  #############################
  
  print("Computing thetas")
  res_winning <- multiomeFate:::peak_mixture_modeling(
    bin_limits = bin_limits,
    bin_midpoints = bin_midpoints, 
    cutmat = cutmat_winning, 
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    bool_freeze_prior = T,
    verbose = 0
  )
  
  res_dying <- multiomeFate:::peak_mixture_modeling(
    bin_limits = bin_limits,
    bin_midpoints = bin_midpoints, 
    cutmat = cutmat_dying, 
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    bool_freeze_prior = F,
    verbose = 0
  )
  
  res_all <- multiomeFate:::peak_mixture_modeling(
    bin_limits = bin_limits,
    bin_midpoints = bin_midpoints, 
    cutmat = cutmat_all, 
    peak_locations = peak_locations,
    peak_prior = peak_prior,
    bool_freeze_prior = F,
    verbose = 0
  )
  
  test_res <- multiomeFate:::peak_testing(
    res_all = res_all,
    res_winning = res_winning,
    res_dying = res_dying
  )
  lrt <- test_res$lrt
  pval <- test_res$pval
  print(paste0("Gene: ", gene, ", -Log10pval = ", round(-log10(pval),2)))
  
  result_list[[i]] <- list(
    num_cells_winning = nrow(cutmat_winning),
    num_cells_dying = nrow(cutmat_dying),
    num_frag_winning = res_winning$num_frags,
    num_frag_dying = res_dying$num_frags,
    theta_winning = res_winning$theta_vec,
    theta_dying = res_dying$theta_vec,
    theta_all = res_all$theta_vec,
    lrt = lrt,
    pval = pval
  )
  
  save(result_list, 
       file = paste0("../../../../out/kevin/Writeup6e/Writeup6e_", treatment, "_peak-entropy_tmp.RData"))
}

save(result_list, treatment,
     date_of_run, session_info,
     winning_cells, dying_cells,
     file = paste0("../../../../out/kevin/Writeup6e/Writeup6e_", treatment, "_peak-entropy.RData"))