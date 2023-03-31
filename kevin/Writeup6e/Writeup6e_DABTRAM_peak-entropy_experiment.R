rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("../Writeup6d/coverage_extractor_singlecell.R")
source("peak_entropy.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "CIS"
source("../Writeup6b/gene_list.R")
source("gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

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

cutmat_winning <- .compute_cutmatrix(
  object = all_data,
  gene = gene,
  cells = winning_cells
)
cutmat_dying <- .compute_cutmatrix(
  object = all_data,
  gene = gene,
  cells = dying_cells
)
cutmat_all <- .compute_cutmatrix(
  object = all_data,
  gene = gene,
  cells = c(winning_cells,dying_cells)
)

peak_mat <- extract_peaks(
  object = all_data,
  gene = gene
)
bin_midpoints <- compute_bin_midpoints(peak_mat)
peak_locations <- compute_peak_locations(peak_mat)
peak_prior <- compute_peak_prior(mat = cutmat_all,
                                 peak_mat = peak_mat)

#############################

res_winning <- peak_mixture_modeling(
  bin_midpoints = bin_midpoints, 
  mat = cutmat_winning, 
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  verbose = F
)
round(res_winning$theta_vec,2)

res_dying <- peak_mixture_modeling(
  bin_midpoints = bin_midpoints, 
  mat = cutmat_dying, 
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  verbose = F
)
round(res_dying$theta_vec,2)

res_all <- peak_mixture_modeling(
  bin_midpoints = bin_midpoints, 
  mat = cutmat_all, 
  peak_locations = peak_locations,
  peak_prior = peak_prior,
  verbose = F
)
round(res_all$theta_vec,2)

lrt <- -2 * (res_all$likelihood_vec[res_all$iter] - (res_winning$likelihood_vec[res_winning$iter] + res_dying$likelihood_vec[res_dying$iter]))
1-stats::pchisq(lrt, df = length(res_dying$theta_vec)-1)

