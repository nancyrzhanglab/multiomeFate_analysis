rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")
source("coverage_extractor_singlecell.R")
source("coverage_extractor_singlecell-plotter.R")

source("../Writeup6b/gene_list.R")
source("gene_list_csc.R")
gene_vec <- sort(unique(c(unlist(keygenes), keygenes_csc)))

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
surviving_lineages <- rownames(tab_mat)[which(tab_mat[,"day10_COCL2"] >= 20)]
winning_idx <- intersect(
  intersect(which(all_data$assigned_lineage %in% surviving_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
dying_idx <- intersect(
  intersect(which(!all_data$assigned_lineage %in% surviving_lineages),
            which(all_data$assigned_posterior >= 0.5)),
  which(all_data$dataset == "day0")
)
length(winning_idx); length(winning_idx)

# initializing
set.seed(10)
rng_idx <- sample(1:length(dying_idx), length(winning_idx))
winning_tracks_list <- vector("list", length = length(gene_vec))
dying_tracks_list <- vector("list", length = length(gene_vec))
peaks_list <- vector("list", length = length(gene_vec))

for(i in 1:length(gene_vec)){
  print(paste0("Gene ", i, " out of ", length(gene_vec), ": ", gene_vec[i]))
  
  window_size <- compute_windowsize(
    object = all_data,
    gene = gene_vec[i],
    base_length = 1000
  )
  if(is.null(window_size)){
    winning_tracks_list[[i]] <- NULL
    dying_tracks_list[[i]] <- NULL
    peaks_list[[i]] <- NULL
    
  } else {
    winning_tracks_list[[i]] <- coverage_extractor_singlecell(
      object = all_data,
      gene = gene_vec[i],
      cells = colnames(all_data)[winning_idx],
      window = window_size
    )
    
    dying_tracks_list[[i]] <- coverage_extractor_singlecell(
      object = all_data,
      gene = gene_vec[i],
      cells = colnames(all_data)[dying_idx[rng_idx]],
      window = window_size
    )
    
    peaks_list[[i]] <- extract_peaks(
      object = all_data,
      gene = gene_vec[i]
    )
  }
}
names(winning_tracks_list) <- gene_vec
names(dying_tracks_list) <- gene_vec
names(peaks_list) <- gene_vec

# remove NULLs
tmp_idx <- which(sapply(peaks_list, function(j){all(is.null(j))}))
if(length(tmp_idx) > 0){
  winning_tracks_list <- winning_tracks_list[-tmp_idx]
  dying_tracks_list <- dying_tracks_list[-tmp_idx]
  peaks_list <- peaks_list[-tmp_idx]
}

###################

for(i in 1:length(peaks_list)){
  print(paste0("Gene ", i, " out of ", length(peaks_list), ": ", names(peaks_list)[i]))
  
  png(paste0("../../../../out/figures/Writeup6d/Writeup6d_COCL2_coverage_", names(peaks_list)[i], ".png"),
      height = 2500, width = 5000, units = "px", res = 300)
  par(mar = c(4,4,4,0.5), mfrow = c(1,2))
  
  plot_coveragetracks(
    cutmat = winning_tracks_list[[i]],
    peaks = peaks_list[[i]],
    max_height = max(winning_tracks_list[[i]]@x)/4,
    main = paste0("COCL2: Winning for ", names(winning_tracks_list)[i]),
  )
  
  plot_coveragetracks(
    cutmat = dying_tracks_list[[i]],
    peaks = peaks_list[[i]],
    max_height = max(dying_tracks_list[[i]]@x)/4,
    main = paste0("COCL2: Losing for ", names(dying_tracks_list)[i], " (Subsampled)"),
  )
  
  graphics.off()
}

