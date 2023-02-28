rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6c/Writeup6c_peak_counter.RData")

countmat_nopeak <- Matrix::Matrix(countmat_nopeak, sparse = T)
quantile(countmat_nopeak@x)

lin_mat <- table(all_data$assigned_lineage, all_data$dataset)
DABTRAM_lineage <- rownames(lin_mat)[which(lin_mat[,"day10_DABTRAM"] >= 20)]
day0_survive_idx <- intersect(which(all_data$assigned_lineage %in% DABTRAM_lineage),
                              which(all_data$dataset == "day0"))

day0_die_idx <- intersect(which(!all_data$assigned_lineage %in% DABTRAM_lineage),
                          which(all_data$dataset == "day0"))

p <- ncol(countmat_nopeak)
pval_vec <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  vec_1 <- as.numeric(countmat_nopeak[day0_survive_idx,j])
  vec_2 <- as.numeric(countmat_nopeak[day0_die_idx,j])
  if(diff(range(vec_1)) == 0 & diff(range(vec_2)) == 0) return(1)
  
  set.seed(10)
  tmp <- stats::wilcox.test(x = vec_1,
                            y = vec_2)
  tmp$p.value
})
names(pval_vec) <- colnames(countmat_nopeak)
# nan_idx <- which(is.nan(pval_vec))

logpval_vec <- -log10(pval_vec)
sort(names(logpval_vec)[which(logpval_vec >= 3)])

source("../Writeup6b/gene_list.R")
gene_vec <- sort(unique(unlist(keygenes)))
logpval_vec[intersect(gene_vec, names(logpval_vec))]

length(which(logpval_vec >= 1))/length(logpval_vec)
length(which(logpval_vec >= 2))/length(logpval_vec)

###################################

dim(countmat_nopeak)
quantile(as.numeric(countmat_nopeak[day0_survive_idx,]))
quantile(as.numeric(countmat_nopeak[day0_die_idx,]))

tmp <- countmat_nopeak[day0_survive_idx,]
quantile(tmp@x, probs = seq(0,1,length.out=11)); length(tmp@x)/prod(dim(tmp))

tmp <- countmat_nopeak[day0_die_idx,]
quantile(tmp@x, probs = seq(0,1,length.out=11)); length(tmp@x)/prod(dim(tmp))

###################################

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

countmat_nopeak_t <- Matrix::t(countmat_nopeak)
nonzero_survive <- sapply(day0_survive_idx, function(i){
  length(.nonzero_col(countmat_nopeak_t, i, bool_value = F))
})
nonzero_die <- sapply(day0_die_idx, function(i){
  length(.nonzero_col(countmat_nopeak_t, i, bool_value = F))
})

quantile(nonzero_survive, probs = seq(0,1,length.out=11))
quantile(nonzero_die, probs = seq(0,1,length.out=11))

###################################

mat <- all_data[["ATAC"]]@counts
mat_t <- Matrix::t(mat)
nonzero_survive2 <- sapply(day0_survive_idx, function(i){
  length(.nonzero_col(mat_t, i, bool_value = F))
})
nonzero_die2 <- sapply(day0_die_idx, function(i){
  length(.nonzero_col(mat_t, i, bool_value = F))
})

quantile(nonzero_survive2, probs = seq(0,1,length.out=11))
quantile(nonzero_die2, probs = seq(0,1,length.out=11))


