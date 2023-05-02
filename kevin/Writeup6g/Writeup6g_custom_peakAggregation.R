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

gene_vec_all <- sort(intersect(gene_vec_all, all_data[["Saver"]]@var.features))
len <- length(gene_vec_all)
result_list <- vector("list", length = len)
names(result_list) <- gene_vec_all

assay <- "ATAC"
extend.downstream <- 5000
extend.upstream <- 5000
sep <- c("-", "-")

for(i in 1:len){
  gene <- gene_vec_all[i]
  print(paste0("Gene ", gene, " (", i, " out of ", len, ")"))
  
  # extract the cutmat
  cutmat <- multiomeFate:::extract_cutmatrix(
    object = all_data,
    gene = gene
  )
  if(all(is.null(cutmat))) next()
  
  # make sure gene exists
  tmp <- Signac::LookupGeneCoords(
    object = all_data,
    gene = gene,
    assay = assay
  )
  if(all(is.null(tmp))) next()
  
  # find the coordinates
  region <- Signac:::FindRegion(
    object = all_data,
    region = gene,
    sep = sep,
    assay = assay,
    extend.upstream = extend.upstream,
    extend.downstream = extend.downstream
  )
  
  # figure out the overlap
  all_atac_peak <- all_data[["ATAC"]]@ranges
  overlap_res <- GenomicRanges::findOverlaps(
    query = all_atac_peak,
    subject = region
  )
  
  if(length(overlap_res) > 0){
    # compute the overlapping regions
    region_gene_peaks <- all_atac_peak[overlap_res@from]
    for(kk in 1:length(region_gene_peaks)){
      region_gene_peaks[kk] <- GenomicRanges::intersect(x = region_gene_peaks[kk],
                                                        y = region)
    }
    
    ranges_obj <- region_gene_peaks@ranges
    if(any(ranges_obj@width <= 150)){
      ranges_obj <- ranges_obj[-which(ranges_obj@width <= 150)]
    }
    if(length(ranges_obj) == 0) next()
    tmp <- cbind(ranges_obj@start, ranges_obj@start + ranges_obj@width - 1)
    colnames(tmp) <- c("start", "end")
    
    colname_vec <- as.numeric(colnames(cutmat))
    chrAct_mat <- sapply(1:nrow(tmp), function(i){
      idx <- intersect(which(colname_vec >= tmp[i,"start"]),
                       which(colname_vec <= tmp[i,"end"]))
      Matrix::rowSums(cutmat[,idx,drop = F])
    })
    
    chromosome_vec <- as.character(region_gene_peaks@seqnames)
    colnames(chrAct_mat) <- sapply(1:nrow(tmp), function(kk){
      paste0(gene, ":", chromosome_vec[kk], ":", tmp[kk,"start"], "-", tmp[kk,"end"])
    })
    
    result_list[[i]] <- Matrix::Matrix(chrAct_mat, sparse = T)
  } 
  
  if(i %% 100 == 0) {
    save(result_list, date_of_run, session_info,
         file = "../../../../out/kevin/Writeup6g/custom_peakAggregation_tmp.RData")
  }
}

save(result_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6g/custom_peakAggregation_tmp.RData")