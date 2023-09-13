rm(list=ls())
library(Seurat); library(Signac)
library(GenomicRanges); library(GenomeInfoDb); library(IRanges)
library(JASPAR2020); library(TFBSTools); library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

treatment <- "COCL2"
Seurat::DefaultAssay(all_data) <- "ATAC"

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

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

keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset %in% c("day0", paste0("week5_", treatment)))] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

motif_focus <- "FOS"

# construct the meta.features for the peaks, which be used for matching
Seurat::DefaultAssay(all_data) <- "ATAC"
mf <- all_data[["ATAC"]]@meta.features[,c("GC.percent", "sequence.length")]
fragment.count <- Matrix::rowSums(all_data[["ATAC"]]@counts)
mf <- cbind(mf, fragment.count)
mf <- scale(mf)
mf <- as.data.frame(mf)
rownames(mf) <- 1:nrow(mf)

##################

motif_matrix <- Signac::GetMotifData(object = all_data, slot = "data")
motif_name_vec <- unlist(all_data[["ATAC"]]@motifs@motif.names)
colnames(motif_matrix) <- motif_name_vec

keep_idx <- grep(motif_focus, motif_name_vec)
motif_matrix <- motif_matrix[,keep_idx]

peak_matrix <- Seurat::GetAssayData(object = all_data, slot = "counts")
niterations <- 100

chromvar_mat <- matrix(NA, nrow = ncol(all_data), ncol = ncol(motif_matrix))
rownames(chromvar_mat) <- colnames(all_data)
colnames(chromvar_mat) <- colnames(motif_matrix)

for(motif_idx in 1:ncol(motif_matrix)){
  print(paste0("Working on motif number ", motif_idx, " out of ", ncol(motif_matrix)))
  motif_name <- colnames(motif_matrix)[motif_idx]
  peak_set <- which(motif_matrix[,motif_name])
  tf_count <- length(peak_set)
  
  # find the background set of peaks
  background_idx <- matrix(NA, nrow = length(peak_set), ncol = niterations)
  colnames(background_idx) <- paste0("Draw:", seq_len(niterations))
  # NOTE: The rows of background_idx are NOT matched to peak_set (i.e., the first row of background_idx
  # are NOT necessarily the matches to peak_set[1])
  
  for(j in 1:niterations){
    if(j %% floor(niterations/10) == 0) cat('*')
    
    set.seed(10*(motif_idx*j))
    tmp <- Signac::MatchRegionStats(
      meta.feature = mf[-peak_set,,drop=F],
      query.feature = mf[peak_set,,drop=F],
      features.match = colnames(mf),
      n = tf_count,
      verbose = F
    )
    if(length(tmp) < tf_count) {
      warning(paste0("Number of requested features in motif number ", motif_idx, 
                     " was too large."))
      tmp <- sample(tmp, size = tf_count, replace = T)
    }
    background_idx[,j] <- as.numeric(tmp)
  }
  
  tf_vec <- Matrix::sparseMatrix(j = peak_set,
                                 i = rep(1, tf_count),
                                 x = 1,
                                 dims = c(1,
                                          nrow(peak_matrix)))
  observed <- as.vector(tf_vec %*% peak_matrix)
  
  niterations <- ncol(background_idx)
  sample_mat <- Matrix::sparseMatrix(j = as.vector(background_idx),
                                     i = rep(seq_len(niterations), each = tf_count),
                                     x = 1,
                                     dims = c(niterations, nrow(peak_matrix)))
  sampled <- as.matrix(sample_mat %*% peak_matrix)
  sampled_mean <- Matrix::colMeans(sampled)
  sampled_sd <- apply(sampled, 2, stats::sd)
  
  z_score <- sapply(1:length(observed), function(i){
    (observed[i] - sampled_mean[i])/sampled_sd[i]
  })
  
  chromvar_mat[names(z_score),motif_name] <- z_score
  
  print("Finished, about to store")
  print(motif_name)
  save(chromvar_mat, motif_focus,
       date_of_run, session_info,
       file = paste0("../../../../out/kevin/Writeup6m/Writeup6m_day0-vs-COCL2week5_chromvar2_", motif_focus, ".RData"))
}

print("Done! :)")