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

gene_vec_all <- sort(intersect(gene_vec_all, rownames(all_data[["RNA"]])))
len <- length(gene_vec_all)
matching_list <- vector("list", length = len)
names(matching_list) <- gene_vec_all

assay <- "ATAC"
extend.downstream <- 5000
extend.upstream <- 5000
sep <- c("-", "-")

for(i in 1:len){
  gene <- gene_vec_all[i]
  print(paste0("Gene ", gene, " (", i, " out of ", len, ")"))
  
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
    overlap_idx <- overlap_res@from
    
    # compute the overlapping regions
    region_gene_peaks <- all_atac_peak[overlap_idx]
    for(kk in 1:length(region_gene_peaks)){
      region_gene_peaks[kk] <- GenomicRanges::intersect(x = region_gene_peaks[kk],
                                                        y = region)
    }
    
    ranges_obj <- region_gene_peaks@ranges
    if(any(ranges_obj@width <= 150)){
      overlap_idx <- overlap_idx[-which(ranges_obj@width <= 150)]
      region_gene_peaks <- region_gene_peaks[-which(ranges_obj@width <= 150)]
    }
    if(length(overlap_idx) == 0) next()
    
    matching_list[[i]] <- list(gene_region = region,
                               overlap_idx = overlap_idx,
                               peak_names = rownames(all_data[["ATAC"]]@counts)[overlap_idx],
                               peak_regions = region_gene_peaks)
  } 
  
  if(i %% 100 == 0) {
    save(matching_list, date_of_run, session_info,
         file = "../../../../out/kevin/Writeup6h/peak-gene-matching.RData")
  }
}

save(matching_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6h/peak-gene-matching.RData")