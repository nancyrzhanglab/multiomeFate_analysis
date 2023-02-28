rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("peak_counter.R")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

gene_vec <- all_data[["Saver"]]@var.features 

extend_width <- 5000
countmat_nopeak <- sapply(1:length(gene_vec), function(i){
  print(paste0("Gene ",i, " out of ", length(gene_vec)))
  
  peak_counter(object = all_data,
               gene = gene_vec[i],
               extend.downstream = extend_width,
               extend.upstream = extend_width)
})
colnames(countmat_nopeak) <- gene_vec

save(countmat_nopeak, extend_width,
     date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6c/Writeup6c_peak_counter_5000.RData")
