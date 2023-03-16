rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6c/Writeup6c_peak_counter_5000_onpeak.RData")
countmat_peak <- countmat_nopeak
load("../../../../out/kevin/Writeup6c/Writeup6c_peak_counter_5000.RData")

countmat_peak <- Matrix::Matrix(countmat_peak, sparse = T)
countmat_nopeak <- Matrix::Matrix(countmat_nopeak, sparse = T)

dim(countmat_peak)
dim(countmat_nopeak)
length(countmat_peak@x)/prod(dim(countmat_peak))
length(countmat_nopeak@x)/prod(dim(countmat_nopeak))

quantile(countmat_peak@x, probs = seq(0,1,length.out=11))
quantile(countmat_nopeak@x, probs = seq(0,1,length.out=11))

####################################

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

.mult_mat_vec <- function(mat, vec){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == ncol(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    mat %*% Matrix::Diagonal(x = vec)
  } else {
    mat * rep(vec, rep(nrow(mat), length(vec)))
  }
}

.mult_vec_mat <- function(vec, mat){
  stopifnot(inherits(mat, c("matrix", "dgCMatrix")), 
            !is.matrix(vec), length(vec) == nrow(mat))
  
  if(inherits(mat, "dgCMatrix")) {
    Matrix::Diagonal(x = vec) %*% mat
  } else {
    vec * mat
  }
}

####################################

p <- ncol(countmat_peak); n <- nrow(countmat_peak)
round(quantile(sapply(1:p, function(j){length(.nonzero_col(countmat_peak, j, F))}), 
               probs = seq(0,1,length.out=11)))
round(quantile(sapply(1:p, function(j){length(.nonzero_col(countmat_nopeak, j, F))}), 
               probs = seq(0,1,length.out=11)))

####################################################

all_data[["offpeak"]] <- Seurat::CreateAssayObject(counts = Matrix::t(countmat_nopeak))

lin_mat <- table(all_data$assigned_lineage, all_data$dataset)
COCL2_lineage <- rownames(lin_mat)[which(lin_mat[,"day10_COCL2"] >= 20)]
day0_survive_idx <- intersect(which(all_data$assigned_lineage %in% COCL2_lineage),
                              which(all_data$dataset == "day0"))
day0_die_idx <- setdiff(intersect(which(!is.na(all_data$assigned_lineage)),
                                  which(all_data$dataset == "day0")),
                        day0_survive_idx)
day0_unknown_idx <- setdiff(which(all_data$dataset == "day0"),
                            c(day0_survive_idx, day0_die_idx))

tmp_vec <- rep(NA, ncol(all_data))
tmp_vec[day0_survive_idx] <- "day0_survive"
tmp_vec[day0_die_idx] <- "day0_die"
tmp_vec[day0_unknown_idx] <- "day0_unknown"
all_data$day0_COCL2_status <- tmp_vec

cell_idx <- sort(unique(c(day0_survive_idx, day0_die_idx, day0_unknown_idx)))
cell_idx2 <- sort(unique(c(day0_survive_idx, day0_die_idx)))
lineage_vec <- all_data$assigned_lineage[cell_idx2]
fitness_vec_onlyday0 <- sapply(1:length(cell_idx2), function(i){
  log1p(lin_mat[lineage_vec[i],"day10_COCL2"])
})
fitness_vec <- rep(NA, ncol(all_data))
fitness_vec[cell_idx2] <- fitness_vec_onlyday0
all_data$day0_COCL2_fitness <- fitness_vec

n <- ncol(all_data)

cell_idx <- sort(unique(c(day0_survive_idx, day0_die_idx, day0_unknown_idx)))
survive_idx <- which(cell_idx %in% day0_survive_idx)
die_idx <- which(cell_idx %in% day0_die_idx)
unknown_idx <- which(cell_idx %in% day0_unknown_idx)

fitness_vec2 <- all_data$day0_COCL2_fitness[cell_idx2]
survive_idx2 <- which(cell_idx2 %in% day0_survive_idx)
die_idx2 <- which(cell_idx2 %in% day0_die_idx)

length(survive_idx) == length(survive_idx2)
length(die_idx) == length(die_idx2)

####################################################

ratio_mat <- (countmat_nopeak)/(countmat_peak+1)
value_vec <- sapply(1:p, function(j){
  num_count <- length(.nonzero_col(countmat_nopeak, j, F))
  log(n/num_count)
})

norm_mat <- .mult_mat_vec(ratio_mat, value_vec)
length(norm_mat@x)/prod(dim(norm_mat))
norm_mat@x <- log1p(norm_mat@x)

norm_mat_t <- Matrix::t(norm_mat)
active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(norm_mat_t, col_idx = i, bool_value = T)
  )
})
peak_plotting_func(active_gene_vec = active_gene_vec, 
                   fitness_vec2 = fitness_vec2, 
                   survive_idx = survive_idx, die_idx = die_idx, unknown_idx = unknown_idx,
                   survive_idx2 = survive_idx2, die_idx2 = die_idx2,
                   filename = "../../../../out/figures/Writeup6c/Writeup6c_COCL2_peakcounter_plotter_custnorm1-offpeak.png",
                   condition_name = "COCL2",
                   main4 = "cust-norm of offpeaks", 
                   xlab = "Custom norm.1 of offpeaks", 
                   ylab4 = "Custom norm.1 of offpeaks")

##############################

countmat_peak_t <- Matrix::t(countmat_peak)

active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(countmat_peak_t, col_idx = i, bool_value = T)
  )
})
peak_plotting_func(active_gene_vec = active_gene_vec, 
                   fitness_vec2 = fitness_vec2, 
                   survive_idx = survive_idx, die_idx = die_idx, unknown_idx = unknown_idx,
                   survive_idx2 = survive_idx2, die_idx2 = die_idx2,
                   filename = "../../../../out/figures/Writeup6c/Writeup6c_COCL2_peakcounter_sum-onpeak.png",
                   condition_name = "COCL2",
                   main4 = "sum of onpeaks", 
                   xlab = "Sum of onpeaks", 
                   ylab4 = "Sum of onpeaks")

##############################

diag_vec <- Matrix::rowSums(countmat_peak)
ratio_mat <- .mult_vec_mat(1/diag_vec, countmat_nopeak)
value_vec <- sapply(1:p, function(j){
  num_count <- length(.nonzero_col(countmat_nopeak, j, T))
  log((n-num_count+1)/(num_count+1))
})
value_vec <- pmax(0, value_vec)
norm_mat <- .mult_mat_vec(ratio_mat, value_vec)*1e4

norm_mat_t <- Matrix::t(norm_mat)
active_gene_vec <- sapply(cell_idx, function(i){
  sum(
    .nonzero_col(norm_mat_t, col_idx = i, bool_value = T)
  )
})
peak_plotting_func(active_gene_vec = active_gene_vec, 
                   fitness_vec2 = fitness_vec2, 
                   survive_idx = survive_idx, die_idx = die_idx, unknown_idx = unknown_idx,
                   survive_idx2 = survive_idx2, die_idx2 = die_idx2,
                   filename = "../../../../out/figures/Writeup6c/Writeup6c_COCL2_peakcounter_plotter_custnorm2-offpeak.png",
                   condition_name = "COCL2",
                   main4 = "cust-norm of offpeaks", 
                   xlab = "Custom norm.2 of offpeaks", 
                   ylab4 = "Custom norm.2 of offpeaks")

