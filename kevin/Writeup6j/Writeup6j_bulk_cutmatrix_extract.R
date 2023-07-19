rm(list=ls())
library(Signac)
library(Rsamtools)
library(GenomicRanges)
library(multiomeFate)
library(rtracklayer)
library(IRanges)

load("../../../../out/kevin/Writeup6g/gene_vec.RData")
source("../Writeup6b/gene_list.R")
source("../Writeup6d/gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

granges_list <- granges_list[gene_vec]
granges_list <- granges_list[sapply(granges_list, length) > 0]
for(i in 1:length(granges_list)){
  if(i %% floor(length(granges_list)/10) == 0) cat('*')
  # from https://github.com/stuart-lab/signac/blob/master/R/utilities.R#L1271
  granges_list[[i]] <- suppressWarnings(Signac::Extend(granges_list[[i]],
                                                       upstream = 5000,
                                                       downstream = 5000))
}

file_prefix <- "~/nzhanglab/data/GoogleDrive_SydneyShafferLab/2019_12_12_EGFR_NGFR_ATAC/Dylans_analysis_scripts/HINT/Merged_Bam_Files/"
file_vec <- c("Mix_Merged_Reads.bam", 
              "EGFR_Merged_Reads.bam", 
              "NGFR_Merged_Reads.bam",
              "EGFR-NGFR-High_Merged_Reads.bam")
names(file_vec) <- c("Mix", "EGFR", "NGFR", "EGFR-NGFR")
pileup_list <- lapply(names(granges_list), function(x){
  lis <- vector("list", length = 4)
  names(lis) <- names(file_vec)
  lis
})
names(pileup_list) <- names(granges_list)
pileupParam <- Rsamtools::PileupParam(max_depth=5000,
                                      distinguish_strands=F, 
                                      distinguish_nucleotides=F)
for(status in names(file_vec)){
  print("=====")
  print(status)
  bamFile <- Rsamtools::BamFile(paste0(file_prefix, file_vec[status]))
  
  for(i in 1:length(granges_list)){
    print(paste0(i, " out of ", length(granges_list)))
    
    gene <- names(granges_list)[i]
    gr <- GenomicRanges::GRanges(seqnames = seqnames(granges_list[[gene]]),
                                 ranges = ranges(granges_list[[gene]]))
    params <- Rsamtools::ScanBamParam(which = gr, what = Rsamtools::scanBamWhat())
    pileup_list[[gene]][[status]] <- Rsamtools::pileup(bamFile, scanBamParam=params,
                                                       pileupParam=pileupParam)
  }
}

# load in the peaks
file_prefix2 <- "~/nzhanglab/data/GoogleDrive_SydneyShafferLab/2019_12_12_EGFR_NGFR_ATAC/Dylans_analysis_scripts/HINT/Merged_Bed_Files/"
file_vec2 <- c("mixed_consensus_peaks.bed", 
               "egfr_consensus_peaks.bed", 
               "ngfr_consensus_peaks.bed",
               "egfr_ngfr_consensus_peaks.bed")
peak_lists <- lapply(file_vec2, function(file){
  import(paste0(file_prefix2, file), format="bed")
})
combined_peaks <- GenomicRanges::reduce(c(peak_lists[[1]], peak_lists[[2]], peak_lists[[3]], peak_lists[[4]]))

# for each gene in names(granges_list), find all the relevant peaks
peak_overlap_list <- lapply(1:length(granges_list), function(i){
  print(paste0(i , " out of ", length(granges_list)))
  
  gene <- names(granges_list)[i]
  overlap_res <- GenomicRanges::findOverlaps(
    query = combined_peaks,
    subject = granges_list[[gene]]
  )
  
  if(length(overlap_res) == 0) return(NULL)
  combined_peaks[overlap_res@from]
})
names(peak_overlap_list) <- names(granges_list)
quantile(sapply(peak_overlap_list, length))
names(peak_overlap_list)[which(sapply(peak_overlap_list, length) == 0)]
