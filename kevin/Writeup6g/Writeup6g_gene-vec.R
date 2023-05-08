rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

gene_vec_all <- rownames(all_data[["RNA"]]@counts)
granges_list <- vector("list", length = length(gene_vec_all))
names(granges_list) <- gene_vec_all

for(i in 1:length(gene_vec_all)){
  if(i %% 100 == 0){
    save(granges_list, date_of_run, session_info,
         file = "../../../../out/kevin/Writeup6g/gene_vec.RData")
  }
  
  print(i)
  gene <- gene_vec_all[i]
  tmp <- Signac::LookupGeneCoords(
    object = all_data,
    gene = gene,
    assay = "ATAC"
  )
  if(!all(is.null(tmp))) granges_list[[i]] <- tmp
}

save(granges_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6g/gene_vec.RData")
