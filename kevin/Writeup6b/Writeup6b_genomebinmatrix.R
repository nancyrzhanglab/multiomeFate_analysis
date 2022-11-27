rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

zz <- all_data[["ATAC"]]
ranges_vec <- zz@ranges@ranges
chrom_vec <- zz@ranges@seqnames
genome <- zz@annotation@seqinfo
seqlengths <- genome@seqlengths
names(seqlengths) <- genome@seqnames

# double check
seqlengths["chr1"]
tmp <- ranges_vec[chrom_vec == "chr1"]
max(tmp@start+tmp@width)

sum(seqlengths)
desired_bin_size <- 50000

genome_bin_matrix <- Signac::GenomeBinMatrix(
  fragments = all_data[["ATAC"]]@fragments,
  genome = seqlengths,
  binsize = desired_bin_size,
  verbose = T)

save(genome_bin_matrix, session_info, date_of_run,
     file = "../../../../out/kevin/Writeup6b/Writeup6b_genomebinmatrix.RData")