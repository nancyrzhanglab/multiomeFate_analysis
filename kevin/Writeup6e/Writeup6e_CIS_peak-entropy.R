rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("../Writeup6d/coverage_extractor_singlecell.R")
source("peak_entropy.R")

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

treatment <- "CIS"

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
  cutmat_all <- rbind(cutmat_winning, cutmat_dying)
  
  peak_mat <- extract_peaks(
    object = all_data,
    gene = gene
  )
  bin_midpoints <- compute_bin_midpoints(peak_mat)
  peak_locations <- compute_peak_locations(peak_mat)
  peak_prior <- compute_peak_prior(mat = cutmat_all,
                                   peak_mat = peak_mat)
  
  #############################
  
  print("Computing thetas")
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
  pval <- 1-stats::pchisq(lrt, df = length(res_dying$theta_vec)-1)
  print(paste0("Gene: ", gene, ", -Log10pval = ", round(-log10(pval),2)))
  
  result_list[[i]] <- list(
    num_cells_winning = nrow(cutmat_winning),
    num_cells_dying = nrow(cutmat_dying),
    num_frag_winning = nrow(res_winning$assignment_mat),
    num_frag_dying = nrow(res_dying$assignment_mat),
    theta_winning = res_winning$theta_vec,
    theta_dying = res_dying$theta_vec,
    theta_all = res_all$theta_vec,
    lrt = lrt,
    pval = pval
  )
}

.compute_entropy <- function(x, tol = 1e-5){
  y <- c(x[1]+x[7], x[2]+x[6], x[3]+x[5], x[4])
  tmp <- sapply(y, function(i){i*log(i)})
  idx <- which(y <= tol)
  if(length(idx)>0) tmp[idx] <- 0
  -sum(tmp)
}
tmp <- sapply(result_list, function(x){
  winning_entropy <- .compute_entropy(x$theta_winning)
  dying_entropy <- .compute_entropy(x$theta_dying)
  c(#x$num_frag_winning, x$num_frag_dying, 
    round(-log10(x$pval+1e-5),2),
    round(winning_entropy - dying_entropy,2))
})
rownames(tmp) <- c(#"frag_winning", "frag_dying", 
  "pval", "entropy_diff")
tmp

round(result_list[["FRZB"]]$theta_winning,2)
round(result_list[["FRZB"]]$theta_dying,2)