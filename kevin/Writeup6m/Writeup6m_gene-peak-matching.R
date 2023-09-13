rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6f/gene_vec.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

gene_vec_all <- sort(intersect(gene_vec_all, rownames(all_data[["RNA"]])))
assay <- "ATAC"
extend.downstream <- 5000
extend.upstream <- 5000

annotations <- Signac::Annotation(object = all_data[["ATAC"]])
tss_positions <- Signac::GetTSSPositions(ranges = annotations)
idx <- which(tss_positions$gene_name %in% gene_vec_all)
tss_positions <- tss_positions[idx]
tss_positions <- tss_positions[order(tss_positions$gene_name)]

gene_vec_all <- tss_positions$gene_name
len <- length(gene_vec_all)
matching_list <- vector("list", length = len)
names(matching_list) <- gene_vec_all

for(i in 1:len){
  gene <- gene_vec_all[i]
  print(paste0("Gene ", gene, " (", i, " out of ", len, ")"))
  
  # find the coordinates
  # FindRegion in https://github.com/stuart-lab/signac/blob/2ad6c3c9c0c8dd31f7e1433b2efd5050d8606f27/R/utilities.R#L1244
  # StringToGRanges in https://github.com/stuart-lab/signac/blob/2ad6c3c9c0c8dd31f7e1433b2efd5050d8606f27/R/utilities.R#L535
  region <- Signac:::FindRegion(
    object = all_data,
    region = tss_positions[i],
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
    region@strand@values <- rep(levels(region@strand)[3], 1)
    region_gene_peaks <- all_atac_peak[overlap_idx]
    for(kk in 1:length(region_gene_peaks)){
      region_gene_peaks[kk] <- GenomicRanges::intersect(x = region_gene_peaks[kk],
                                                        y = region)
    }
    
    ranges_obj <- region_gene_peaks@ranges
    if(any(ranges_obj@width <= 50)){
      overlap_idx <- overlap_idx[-which(ranges_obj@width <= 50)]
      region_gene_peaks <- region_gene_peaks[-which(ranges_obj@width <= 50)]
    }
    if(length(overlap_idx) == 0) next()
    
    matching_list[[i]] <- list(gene_tss = tss_positions[i],
                               overlap_idx = overlap_idx,
                               peak_names = rownames(all_data[["ATAC"]]@counts)[overlap_idx],
                               peak_regions = region_gene_peaks)
  } 
  
  if(i %% 100 == 0) {
    save(matching_list, date_of_run, session_info,
         file = "../../../../out/kevin/Writeup6m/peak-gene-matching.RData")
  }
}

notes <- paste0("Unlike Writeup6h_gene-peak-matching.R, this version: 1) Uses ",
                "TSS locations via Signac::GetTSSPositions instead of ",
                "Signac::LookupGeneCoords, and 2) throws out peaks that are ",
                "less than 50 bp instead of 150 bp.")

save(matching_list, notes, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup6m/peak-gene-matching.RData")