rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6i/day0_binomial_extract-DABTRAM_tmp.RData")
source("fisher_exact_test.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

tmp <- sapply(result_list, function(x){all(!is.null(x))})
result_list <- result_list[which(tmp)]

pvalue_func <- function(mat){
  bin_win <- mat[1,]
  bin_die <- mat[2,]
  
  if(all(bin_win[c("Mid", "Close")] == 0) | all(bin_die[c("Mid", "Close")] == 0)){
    val1 <- NA
  } else {
    count_group1 <- bin_win[c("Close", "Mid")]
    count_group2 <- bin_die[c("Close", "Mid")]
    res1 <- exact_fisher_pvalue(count_group1 = count_group1,
                                count_group2 = count_group2)
    val1 <- res1$log10prob_val
  }
  
  if(all(bin_win[c("Far", "Close")] == 0) | all(bin_die[c("Far", "Close")] == 0)){
    val2 <- NA
  } else {
    count_group1 <- bin_win[c("Close", "Far")]
    count_group2 <- bin_die[c("Close", "Far")]
    res2 <- exact_fisher_pvalue(count_group1 = count_group1,
                                count_group2 = count_group2)
    val2 <- res2$log10prob_val
  }
  
  c(val1, val2)
}

len <- length(result_list)
pvalue_mat <- sapply(1:len, function(i){
  if(i %% floor(len/10) == 0) cat('*')
  if(length(result_list[[i]]) == 0) return(rep(NA,2))
  pvalue_func(result_list[[i]])
})
colnames(pvalue_mat) <- names(result_list)
rownames(pvalue_mat) <- c("Mid", "Far")

vec1 <- pvalue_mat[1,]
vec1 <- vec1[!is.na(vec1)]
vec1 <- 10^(-vec1)

vec2 <- pvalue_mat[2,]
vec2 <- vec2[!is.na(vec2)]
vec2 <- 10^(-vec2)

png(paste0("../../../../out/figures/Writeup6i/Writeup6i_day0_binomial-pvalue_DABTRAM.png"),
    height = 1200, width = 2500, res = 300, units = "px")
par(mfrow = c(1,2))
hist(vec1, breaks = 15, xlab = "P-value", main = "Mid vs. Close (DABTRAM)")

hist(vec2, breaks = 15, xlab = "P-value", main = "Far vs. Close (DABTRAM)")
graphics.off()

vec1b <- stats::p.adjust(vec1, method = "BH")
vec1b[vec1b <= .05]

vec2b <- stats::p.adjust(vec2, method = "BH")
vec2b[vec2b <= .05]
