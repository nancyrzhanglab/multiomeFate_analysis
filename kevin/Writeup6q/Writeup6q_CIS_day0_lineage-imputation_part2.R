rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

treatment <- "CIS"
day_early <- "day0"
day_later <- "day10"
day_early_full <- day_early
day_later_full <- paste0(day_later, "_", treatment)

load(paste0("../../../../out/kevin/Writeup6q/Writeup6q_", treatment, "_", day_early, "_lineage-imputation.RData"))

###################

train_mat <- sapply(loocv_fit_list, function(x){
  x$train_loglik
})
Matrix::rowMeans(train_mat)

test_mat <- sapply(loocv_fit_list, function(x){
  x$test_loglik
})
loglik_mean <- Matrix::rowMeans(test_mat)
loglik_mean
