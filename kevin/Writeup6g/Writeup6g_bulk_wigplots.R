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

pdf("../../../../out/figures/Writeup6g/Writeup6g_bulk_wigplot.pdf", 
    onefile = T, width = 5, height = 8)
par(mfrow = c(4,1), mar = c(4,4,4,0.5))
for(i in 1:length(pileup_list)){
  print(paste0(i, " out of ", length(pileup_list)))
  
  gene <- names(pileup_list)[i]
  x_vec <- (granges_list[[gene]]@ranges@start):(granges_list[[gene]]@ranges@start + granges_list[[gene]]@ranges@width - 1)
  
  for(status in names(file_vec)){
    x_vec2 <- pileup_list[[gene]][[status]][,"pos"]
    y_vec2 <- pileup_list[[gene]][[status]][,"count"]
    y_vec <- rep(0, length(x_vec))
    idx <- x_vec %in% x_vec2
    y_vec[idx] <- y_vec2
    y_vec <- y_vec/max(y_vec)
    
    plot(NA, xlim = range(x_vec), ylim = range(y_vec), 
         xlab = "Basepair", ylab = "Rescaled count",
         main = paste0(gene, ": ", status), bty="n")
    for(y in seq(0,1,length.out=11)){
      lines(x = range(x_vec), y = rep(y, 2), 
            lty = 2, lwd = 0.5, col = "lightgray")
    }
    
    #plot peaks
    peak_mat <- peak_overlap_list[[gene]]
    if(length(peak_mat) > 0){
      for(i in 1:length(peak_mat)){
        tmp <- peak_mat[i]@ranges
        val <- c(tmp@start, tmp@start+tmp@width-1)
        graphics::polygon(
          x = c(val[c(1,2,2,1)]),
          y = c(2,2,-2,-2),
          border = NA,
          density = NULL,
          col = rgb(255, 102, 102, alpha = 0.2*255, maxColorValue = 255)
        )
      }
    }
    
    polygon(x = c(x_vec[1], x_vec, x_vec[length(x_vec)]),
            y = c(0, y_vec, 0),
            col = "gray25",
            border = NA,
            density = NULL)
  }
}
dev.off()