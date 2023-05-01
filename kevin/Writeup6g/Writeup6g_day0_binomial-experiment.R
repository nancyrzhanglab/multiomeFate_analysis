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

#####################################

gene <- "MYL6"
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
peak_prior <-  multiomeFate:::compute_peak_prior(cutmat = cutmat_all,
                                                 peak_mat = peak_mat)
# remove peaks with less than 0.05
if(any(peak_prior <= 0.05)){
  idx <- which(peak_prior <= 0.05)
  print(paste0("Removed " , length(idx), " peaks"))
  peak_locations <- peak_locations[-idx]
  peak_prior <- peak_prior[-idx]
  peak_prior <- peak_prior/sum(peak_prior)
}


###############

# start to formulate this in terms of the binomial distribution
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
names(bin_win) <- c("Close", "Mid", "Far")
bin_die <- sapply(1:3, function(k){
  sum(collapse_die[which(assignment_vec == k)])
})
names(bin_die) <- c("Close", "Mid", "Far")

# let's first try a prop.test
x_mat <- rbind(bin_win[c("Mid", "Close")], bin_die[c("Mid", "Close")])
res1a <- stats::prop.test(x = x_mat, alternative = "greater")
x_mat <- rbind(bin_win[c("Far", "Close")], bin_die[c("Far", "Close")])
res2a <- stats::prop.test(x = x_mat, alternative = "greater")
res1a; res2a

# now let's try a fisher exact test
# this one is a bit more sus since it doesn't account for the variability in the null
fisher_exact_test <- function(obs_success, n, null_p){
  1 - stats::pbinom(obs_success, size = n, prob = null_p)
}
res1b <- fisher_exact_test(obs_success = bin_win["Mid"],
                           n = sum(bin_win[c("Mid", "Close")]),
                           null_p = bin_die["Mid"]/sum(bin_die[c("Mid", "Close")]))
res2b <- fisher_exact_test(obs_success = bin_win["Far"],
                           n = sum(bin_win[c("Far", "Close")]),
                           null_p = bin_die["Far"]/sum(bin_die[c("Far", "Close")]))


