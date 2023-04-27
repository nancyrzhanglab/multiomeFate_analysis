rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

gene_vec_all <- rownames(all_data[["Saver"]]@data)
tmp <- sapply(1:length(gene_vec_all), function(i){
  print(i)
  # if(i %% floor(length(gene_vec_all)/10) == 0) cat('*')
  gene <- gene_vec_all[i]
  tmp <- Signac::LookupGeneCoords(
    object = all_data,
    gene = gene,
    assay = "ATAC"
  )
  if(all(is.null(tmp))) return(F) else T
})
gene_vec_all <- gene_vec_all[which(tmp)]
gene_vec_all <- sort(gene_vec_all)

save(gene_vec_all, date_of_run, session_info,
     file = paste0("../../../../out/kevin/Writeup6f/gene_vec_shorten.RData"))