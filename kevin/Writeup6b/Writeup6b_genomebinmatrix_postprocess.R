rm(list=ls())
library(Seurat)
library(Signac)

load("../../../../out/kevin/Writeup6b/Writeup6b_genomebinmatrix.RData")
load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

dim(genome_bin_matrix)
genome_bin_matrix[1:5,1:5]
sum(genome_bin_matrix@x)

######################
# let's test one random chromosome region
######################
zz <- all_data[["ATAC"]]
ranges_vec <- zz@ranges@ranges
chrom_vec <- zz@ranges@seqnames
length(chrom_vec) == length(ranges_vec)
chrom <- "chr4"
chrom_idx <- grep(chrom, rownames(genome_bin_matrix))
chrom_val <- Matrix::rowSums(genome_bin_matrix[chrom_idx,])
chrom_val <- sort(chrom_val, decreasing = T)
head(chrom_val)

# find a suitable chromosome region
# chrom_region <- names(chrom_val)[round(length(chrom_val)/4)]
chrom_region <- names(chrom_val)[1]
start <- strsplit(chrom_region, "-")[[1]][2]
end <- strsplit(chrom_region, "-")[[1]][3]
chrom_computed <- genome_bin_matrix[chrom_region,]

# find it's relevant rows in ATAC
idx1 <- which(as.character(chrom_vec) == chrom)
ranges_mat <- as.matrix(ranges_vec)
idx2 <- which(ranges_mat[,1] <= start)
idx3 <- which(ranges_mat[,1]+ranges_mat[,2] <= end)
idx4 <- which(ranges_mat[,1] >= start)
idx5 <- which(ranges_mat[,1]+ranges_mat[,2] >= end)

# idx_final <- setdiff(idx1, intersect(idx4, idx5))
# length(idx_final)
# idx_final <- setdiff(idx_final, intersect(idx2, idx3))
# length(idx_final)

idx_final <- intersect(idx1, intersect(idx4, idx3))
length(idx_final)

chrom_target <- Matrix::colSums(all_data[["ATAC"]]@counts[idx_final,])

stats::cor(chrom_computed, chrom_target)
tmp <- chrom_computed - chrom_target
quantile(chrom_computed, probs = seq(0,1,length.out=11))
quantile(chrom_target, probs = seq(0,1,length.out=11))
quantile(tmp, probs = seq(0,1,length.out=11))
round(100*length(which(tmp != 0))/length(tmp))

table(chrom_computed, tmp)




