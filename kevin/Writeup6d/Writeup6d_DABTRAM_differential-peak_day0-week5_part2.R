rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
load("../../../../out/kevin/Writeup6d/Writeup6d_DABTRAM_differential-peak_day0-week5.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

treatment <- "DABTRAM"
Seurat::DefaultAssay(all_data) <- "ATAC"

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_names <- rownames(tab_mat)[which(tab_mat[,paste0("week5_", treatment)] >= 10)]
cell_names1 <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names)]
cell_names2 <- colnames(all_data)[which(all_data$dataset == "day0")]
cell_names_winning <- intersect(cell_names1, cell_names2)
cell_names_losing <- setdiff(cell_names2, cell_names1)
ident_vec <- rep(NA, ncol(all_data))
names(ident_vec) <- colnames(all_data)
ident_vec[cell_names_winning] <- paste0("day0_win_", treatment)
ident_vec[cell_names_losing] <- paste0("day0_lose_", treatment)
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"
table(Seurat::Idents(all_data))

#########################

idx <- intersect(which(de_res[,"p_val"]<=1e-2), 
                 which(de_res[,"avg_log2FC"] > 0))
head(de_res[idx,],10)
length(idx)
peak_names <- sort(rownames(de_res)[idx])
zz <- all_data[["ATAC"]]
peak_idx <- which(rownames(all_data[["ATAC"]]) %in% peak_names)
pos_ranges <- all_data[["ATAC"]]@ranges[peak_idx]
length(pos_ranges)
sort(unique(pos_ranges@seqnames))

idx <- intersect(which(de_res[,"p_val"]<=1e-2), 
                 which(de_res[,"avg_log2FC"] < 0))
head(de_res[idx,],10)
length(idx)
peak_names <- sort(rownames(de_res)[idx])
zz <- all_data[["ATAC"]]
peak_idx <- which(rownames(all_data[["ATAC"]]) %in% peak_names)
neg_ranges <- all_data[["ATAC"]]@ranges[peak_idx]
length(neg_ranges)
sort(unique(neg_ranges@seqnames))

.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, c("dgCMatrix", "lgCMatrix")), col_idx %% 1 == 0,
            col_idx > 0, col_idx <= ncol(mat))
  
  val1 <- mat@p[col_idx]
  val2 <- mat@p[col_idx+1]
  
  if(val1 == val2) return(numeric(0))
  if(bool_value){
    # return the value
    mat@x[(val1+1):val2]
  } else {
    # return the row index
    mat@i[(val1+1):val2]+1
  }
}

cell_idx <- which(all_data$ident %in% c(paste0("day0_win_", treatment), paste0("day0_lose_", treatment)))
tmp_mat <- all_data[["ATAC"]]@counts[,cell_idx]
tmp_mat <- Matrix::t(tmp_mat)
num_cells_per_peak <- sapply(1:ncol(tmp_mat), function(j){
  length(.nonzero_col(tmp_mat, j, F))
})
bg_peak_idx <- which(num_cells_per_peak >= 0.1*length(cell_idx))
bg_ranges <- all_data[["ATAC"]]@ranges[bg_peak_idx]
length(bg_ranges)
sort(unique(bg_ranges@seqnames))

# now to make three tab-separated files, one with the enriched regions (in each direction) and one with all other regions
write_peakfile <- function(range_obj, file){
  tmp <- as.character(range_obj@seqnames)
  rm_idx <- which(!tmp %in% paste0("chr", c(1:22)))
  if(length(rm_idx) > 0){
    range_obj <- range_obj[-rm_idx]
  }
  
  chr_vec <- as.character(range_obj@seqnames)
  start_vec <- range_obj@ranges@start
  end_vec <- range_obj@ranges@start + range_obj@ranges@width - 1
  
  n <- length(range_obj)
  fileConn <- file(file)
  vec <- sapply(1:n, function(i){
    id <- paste0(chr_vec[i], ":", start_vec[i], "-", end_vec[i])
    paste(c(id, chr_vec[i], start_vec[i], end_vec[i], "0"), collapse = "\t")
  })
  writeLines(vec, fileConn, sep = "\n")
  close(fileConn)
  
  invisible()
}

write_peakfile(pos_ranges, 
               file = "../../../../out/kevin/Writeup6d/Writeup6d_DABTRAM_day0-week5_differential_pospeaks.txt")

write_peakfile(neg_ranges, 
               file = "../../../../out/kevin/Writeup6d/Writeup6d_DABTRAM_day0-week5_differential_negpeaks.txt")

write_peakfile(bg_ranges, 
               file = "../../../../out/kevin/Writeup6d/Writeup6d_DABTRAM_day0-week5_differential_bgpeaks.txt")

