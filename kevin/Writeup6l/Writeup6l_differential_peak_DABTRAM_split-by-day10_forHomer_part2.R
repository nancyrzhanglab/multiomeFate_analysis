rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

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

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
keep_vec <- rep(FALSE, ncol(all_data))
keep_vec[which(all_data$dataset == "day0")] <- TRUE
all_data$keep <- keep_vec
all_data <- subset(all_data, keep == TRUE)

################

print("Creating partitions")
treatment <- "DABTRAM"
Seurat::DefaultAssay(all_data) <- "ATAC"

lineage_names_win <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] >= 10)]
cell_names_win <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names_win)]
lineage_names_lose <- rownames(tab_mat)[which(tab_mat[,paste0("day10_", treatment)] == 0)]
cell_names_lose <- colnames(all_data)[which(all_data$assigned_lineage %in% lineage_names_lose)]
ident_vec <- rep(NA, ncol(all_data))
names(ident_vec) <- colnames(all_data)
ident_vec[cell_names_win] <- "day0_win"
ident_vec[cell_names_lose] <- "day0_lose"
all_data$ident <- ident_vec
Seurat::Idents(all_data) <- "ident"
table(Seurat::Idents(all_data))

###############################

print("Adding motifs")
pfm <- getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# add motif information
all_data <- Signac::AddMotifs(
  object = all_data,
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
  pfm = pfm
)

print("Computing peak statistics")
all_data <- Signac::RegionStats(
  object = all_data, 
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
  assay = "ATAC")
head(all_data[["ATAC"]]@meta.features)

###############################

print("Finding differentially-accessible peaks")
# from https://stuartlab.org/signac/articles/pbmc_vignette.html#find-differentially-accessible-peaks-between-clusters
set.seed(10)
de_res <- Seurat::FindMarkers(
  object = all_data,
  ident.1 = "day0_win",
  ident.2 = "day0_lose",
  test.use = 'LR',
  only.pos = FALSE,
  latent.vars = 'nCount_ATAC',
  min.pct = 0.05,
  verbose = T
)

save(date_of_run, session_info,
     de_res,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_differential_peak_DABTRAM_split-by-day10_part2.RData")

#################


print("Partitioning pos/neg peaks")
idx <- intersect(which(de_res[,"p_val"] <= 1e-2), 
                 which(de_res[,"avg_log2FC"] > 0))
print(paste0("Number of positive enriched peaks: ", length(idx)))
pos_names <- sort(rownames(de_res)[idx])
length(pos_names)

idx <- intersect(which(de_res[,"p_val"] <= 1e-2), 
                 which(de_res[,"avg_log2FC"] < 0))
print(paste0("Number of negatively enriched peaks: ", length(idx)))
neg_names <- sort(rownames(de_res)[idx])
length(neg_names)

save(date_of_run, session_info,
     pos_names, neg_names, de_res, 
     file = "../../../../out/kevin/Writeup6l/Writeup6l_differential_peak_DABTRAM_split-by-day10_part2.RData")

###############################

# https://stuartlab.org/signac/reference/findmotifs
# https://stuartlab.org/signac/articles/motif_vignette
print("Grabbing open peaks")
open_peaks <- Signac::AccessiblePeaks(
  all_data, 
  idents = c("day0_win", "day0_lose")
)

meta_feature <- Seurat::GetAssayData(
  all_data, 
  assay = "ATAC", 
  slot = "meta.features")

# remove non-standard chromatin annotations
na_idx <- which(apply(meta_feature, 1, function(x){any(is.na(x))}))
if(length(na_idx) > 0){
  meta_feature <- meta_feature[-na_idx,,drop = F]
}
open_peaks <- open_peaks[open_peaks %in% rownames(meta_feature)]
pos_names <- pos_names[pos_names %in% rownames(meta_feature)]
neg_names <- neg_names[neg_names %in% rownames(meta_feature)]

save(date_of_run, session_info,
     pos_names, neg_names, de_res, 
     open_peaks, meta_feature,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_differential_peak_DABTRAM_split-by-day10_part2.RData")

print("Finding pos matching peaks")
# match the overall GC content in the peak set
peaks_matched <- Signac::MatchRegionStats(
  meta.feature = meta_feature[open_peaks, ],
  query.feature = meta_feature[pos_names, ],
  n = 50000
)
print("Finding pos motifs")
set.seed(10)
pos_motif_de <- Signac::FindMotifs(
  object = all_data,
  features = pos_names,
  background = peaks_matched
)

print("Finding neg matching peaks")
# match the overall GC content in the peak set
meta_feature <- Seurat::GetAssayData(
  all_data, 
  assay = "ATAC", 
  slot = "meta.features")
peaks_matched <- Signac::MatchRegionStats(
  meta.feature = meta_feature[open_peaks, ],
  query.feature = meta_feature[neg_names, ],
  n = 50000
)
print("Finding neg motifs")
set.seed(10)
neg_motif_de <- Signac::FindMotifs(
  object = all_data,
  features = neg_names,
  background = peaks_matched
)

print("Saving")
save(date_of_run, session_info,
     pos_names, neg_names, de_res,
     pos_motif_de, neg_motif_de,
     open_peaks, meta_feature,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_differential_peak_DABTRAM_split-by-day10_part2.RData")

##################

pos_motif_de[intersect(which(pos_motif_de$p.adjust <= 0.05), 1:100),]
neg_motif_de[intersect(which(neg_motif_de$p.adjust <= 0.05), 1:100),]

pos_motif <- sort(pos_motif_de[which(pos_motif_de$p.adjust <= 0.05),"motif.name"])
neg_motif <- sort(neg_motif_de[which(neg_motif_de$p.adjust <= 0.05),"motif.name"])
common_motif <- sort(intersect(pos_motif, neg_motif))
pos_motif <- setdiff(pos_motif, common_motif)
neg_motif <- setdiff(neg_motif, common_motif)

pos_motif
length(pos_motif)

neg_motif
length(neg_motif)

