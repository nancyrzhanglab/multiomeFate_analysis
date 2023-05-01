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

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
treatment_vec <- c("CIS", "COCL2", "DABTRAM")

result_list <- vector("list", length = 3)
names(result_list) <- treatment_vec

for(treatment in treatment_vec){
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
  winning_cells <- colnames(all_data)[winning_idx]
  dying_cells <- colnames(all_data)[dying_idx]
  
  #############
  result_list[[treatment]] <- lapply(1:length(gene_vec_all), function(i){
    print(paste0(i, " of ", length(gene_vec_all), "in ", treatment))
    gene <- gene_vec_all[i]
    
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
    peak_prior <-  multiomeFate:::compute_peak_prior(cutmat = cutmat_all,
                                                     peak_mat = peak_mat)
    # remove peaks with less than 0.05
    if(any(peak_prior <= 0.05)){
      if(length(which(peak_prior >= 0.05)) == 0) return(NULL)
      idx <- which(peak_prior <= 0.05)
      print(paste0("Removed " , length(idx), " peaks"))
      peak_locations <- peak_locations[-idx]
      peak_prior <- peak_prior[-idx]
      peak_prior <- peak_prior/sum(peak_prior)
    }

    collapse_win <- Matrix::colSums(cutmat_winning)
    collapse_die <- Matrix::colSums(cutmat_dying)
    
    # compute a vector of length 3: the number more than 2x from peak-midpoint, .5x-2x peak midpoint, closer than .5x peak midpoint
    assignment_vec <- rep(3, ncol(cutmat_winning))
    names(assignment_vec) <- colnames(cutmat_winning)
    
    colname_vec <- as.numeric(names(assignment_vec))
    for(i in 1:nrow(peak_mat)){
      tmp <- peak_mat[i,]
      start <- tmp[1]; end <- tmp[2]
      midpoint <- c(start+end)/2
      len <- end - start
      # for <=2x
      vec <- midpoint + c(-2,2)*len
      idx <- intersect(which(colname_vec >= vec[1]), which(colname_vec <= vec[2]))
      assignment_vec[idx] <- 2
      
      # for <=.5x
      vec <- midpoint + c(-.5,.5)*len
      idx <- intersect(which(colname_vec >= vec[1]), which(colname_vec <= vec[2]))
      assignment_vec[idx] <- 1
    }
    
    bin_win <- sapply(1:3, function(k){
      sum(collapse_win[which(assignment_vec == k)])
    })
    bin_die <- sapply(1:3, function(k){
      sum(collapse_die[which(assignment_vec == k)])
    })
    
    res <- rbind(bin_win, bin_die)
    colnames(res) <- c("Close", "Mid", "Far")
    rownames(res) <- c(paste0("Win_day0_", treatment),
                       paste0("Die_day0_", treatment))
    res
  })
  save(result_list, date_of_run, session_info,
       file = "../../../../out/kevin/Writeup6g/day0_binomial_extract-all_tmp.RData")
  
  names(result_list[[treatment]]) <- gene_vec_all
  save(result_list, date_of_run, session_info,
       file = "../../../../out/kevin/Writeup6g/day0_binomial_extract-all_tmp.RData")
}