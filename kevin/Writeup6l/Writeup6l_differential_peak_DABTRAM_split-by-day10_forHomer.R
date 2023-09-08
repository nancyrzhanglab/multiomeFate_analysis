rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

print("Simplifying dataset")
Seurat::DefaultAssay(all_data) <- "Saver"
all_data[["geneActivity"]] <- NULL

all_data[["pca"]] <- NULL
all_data[["umap"]] <- NULL
all_data[["lsi"]] <- NULL
all_data[["atac.umap"]] <- NULL
all_data[["saverpca"]] <- NULL
all_data[["fasttopic_CIS"]] <- NULL
all_data[["fasttopic_COCL2"]] <- NULL
all_data[["fasttopic_DABTRAM"]] <- NULL
all_data[["common_tcca"]] <- NULL
all_data[["distinct1_tcca"]] <- NULL
all_data[["distinct2_tcca"]] <- NULL

print("Removing cells without a barcode")
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[intersect(which(!is.na(all_data$assigned_lineage)),
                   which(all_data$assigned_posterior >= 0.5))] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

################

treatment <- "DABTRAM"
Seurat::DefaultAssay(all_data) <- "ATAC"

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
lineage_names_win <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 20)]
cell_names_win <- colnames(all_data)[intersect(
  which(all_data$assigned_lineage %in% lineage_names_win),
  which(all_data$dataset == "day0")
)]
lineage_names_lose <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] == 0)]
cell_names_lose <- colnames(all_data)[intersect(
  which(all_data$assigned_lineage %in% lineage_names_lose),
  which(all_data$dataset == "day0")
)]
ident_vec <- rep(NA, ncol(all_data))
names(ident_vec) <- colnames(all_data)
ident_vec[cell_names_win] <- "day0_win"
ident_vec[cell_names_lose] <- "day0_lose"
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"
table(Seurat::Idents(all_data))

# from https://stuartlab.org/signac/articles/pbmc_vignette.html#find-differentially-accessible-peaks-between-clusters
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = paste0("day0_win_", treatment),
  ident.2 = paste0("day0_lose_", treatment),
  test.use = 'LR',
  latent.vars = 'nCount_ATAC',
  verbose = T
)

save(date_of_run, session_info, de_res,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_DABTRAM-split-by-day10_differential-peak.RData")

#################

idx <- intersect(which(de_res[,"p_val"] <= 1e-2), 
                 which(de_res[,"avg_log2FC"] > 0))
print(paste0("Number of positive enriched peaks: ", length(idx)))
peak_names <- sort(rownames(de_res)[idx])
peak_idx <- which(rownames(all_data[["ATAC"]]) %in% peak_names)
pos_ranges <- all_data[["ATAC"]]@ranges[peak_idx]

idx <- intersect(which(de_res[,"p_val"] <= 1e-2), 
                 which(de_res[,"avg_log2FC"] < 0))
print(paste0("Number of negatively enriched peaks: ", length(idx)))
peak_names <- sort(rownames(de_res)[idx])
peak_idx <- which(rownames(all_data[["ATAC"]]) %in% peak_names)
neg_ranges <- all_data[["ATAC"]]@ranges[peak_idx]

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
bg_peak_idx <- which(num_cells_per_peak >= 0.05*length(cell_idx))
bg_ranges <- all_data[["ATAC"]]@ranges[bg_peak_idx]

# now to make three tab-separated files, one with the enriched regions (in each direction) and one with all other regions
write_peakfile <- function(range_obj, file){
  tmp <- as.character(range_obj@seqnames)
  rm_idx <- which(tmp %in% c("chrX", "GL000195.1", "GL000205.2", "GL000219.1", "KI270713.1"))
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
               file = "../../../../out/kevin/Writeup6l/Writeup6l_DABTRAM-split-by-day10_differential_pospeaks.txt")

write_peakfile(neg_ranges, 
               file = "../../../../out/kevin/Writeup6l/Writeup6l_DABTRAM-split-by-day10_differential_negpeaks.txt")

write_peakfile(bg_ranges, 
               file = "../../../../out/kevin/Writeup6l/Writeup6l_DABTRAM-split-by-day10_differential_bgpeaks.txt")

