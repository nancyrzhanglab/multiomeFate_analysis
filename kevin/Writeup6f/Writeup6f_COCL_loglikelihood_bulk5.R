rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6f/gene_vec_shorten.RData")
gene_vec_all <- sort(gene_vec_all)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
bulk_number <- 5
bulk_total <- 5

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

gene_len <- length(gene_vec_all)
num_per_bulk <- ceiling(gene_len/bulk_total)
gene_vec <- gene_vec_all[(num_per_bulk*(bulk_number-1)+1):(min(gene_len, num_per_bulk*bulk_number))]
result_list <- vector("list", length = length(gene_vec))
names(result_list) <- gene_vec

###############################

for(i in 1:length(gene_vec)){
  gene <- gene_vec[i]
  print(gene)
  
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
  
  if(length(cutmat_winning@x) <= 5 | length(cutmat_dying@x) <= 5) next()
  
  peak_mat <- multiomeFate:::extract_peaks(
    object = all_data,
    gene = gene
  )
  peak_locations <- multiomeFate:::compute_peak_locations(peak_mat)
  peak_prior <-  multiomeFate:::compute_peak_prior(cutmat = cutmat_all,
                                                   peak_mat = peak_mat)
  peak_width <- stats::median(apply(peak_mat, 1, diff))
  
  # remove peaks with less than a prior of 0.01
  if(any(peak_prior <= 0.05)){
    idx <- which(peak_prior <= 0.05)
    peak_locations <- peak_locations[-idx]
    peak_prior <- peak_prior[-idx]
    peak_prior <- peak_prior/sum(peak_prior)
  }
  
  set.seed(i)
  res_win <- multiomeFate:::peak_mixture_modeling(cutmat = cutmat_dying,
                                                  peak_locations = peak_locations,
                                                  peak_prior = peak_prior,
                                                  peak_width = peak_width,
                                                  return_dist_mat = F,
                                                  max_iter = 100,
                                                  verbose = 0)
  res_die <- multiomeFate:::peak_mixture_modeling(cutmat = cutmat_winning,
                                                  peak_locations = peak_locations,
                                                  peak_prior = peak_prior,
                                                  peak_width = peak_width,
                                                  return_dist_mat = F,
                                                  max_iter = 100,
                                                  verbose = 0)
  res_both <- multiomeFate:::peak_mixture_modeling(cutmat = rbind(cutmat_winning,
                                                                  cutmat_dying),
                                                   peak_locations = peak_locations,
                                                   peak_prior = peak_prior,
                                                   peak_width = peak_width,
                                                   return_dist_mat = F,
                                                   max_iter = 100,
                                                   verbose = 0)
  
  result_list[[i]] <- list(res_win = res_win,
                           res_die = res_die,
                           res_both = res_both)
  
  save(result_list, date_of_run, session_info,
       file = paste0("../../../../out/kevin/Writeup6f/COCL_loglikelihood_bulk", bulk_number, "of", bulk_total, "_tmp.RData"))
}