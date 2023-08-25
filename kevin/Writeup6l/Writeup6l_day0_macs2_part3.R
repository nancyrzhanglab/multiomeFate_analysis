rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

print("Loading data")
load("../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2.RData") #load: peaks
load("../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2_part2.RData") #load: feature_mat
load("../../../../out/kevin/Writeup6l/Writeup6l_day0_chromvar_lightweight_noATAC.RData") #load: all_data (with no ATAC)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

print("Clearing elements")
all_data[["chromvar"]] <- NULL
all_data$nCount_ATAC <- NULL
all_data$nFeature_ATAC <- NULL
all_data$nucleosome_signal <- NULL
all_data$nucleosome_percentile <- NULL
all_data$TSS.enrichment <- NULL
all_data$TSS.percentile <- NULL
all_data$high.tss <- NULL
all_data$blacklist_fraction <- NULL
all_data$passed_filters <- NULL
all_data$peak_region_fragments <- NULL
all_data$pct_reads_in_peaks <- NULL
all_data$nCount_geneActivity <- NULL
all_data$nFeature_geneActivity <- NULL
all_data$nCount_spliced <- NULL
all_data$nFeature_spliced <- NULL
all_data$nCount_unspliced <- NULL
all_data$nFeature_unspliced <- NULL

print("Tidying up feature and fragments")
colnames(feature_mat) <- paste0("day0_", colnames(feature_mat))
all(colnames(feature_mat) %in% colnames(all_data))

cells <- colnames(all_data)
cells <- sapply(strsplit(cells, split = "_"), function(x){x[2]})
fragments <- Signac::CreateFragmentObject("~/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time0/outs/atac_fragments.tsv.gz",
                                          cells = cells)

print("Setting up annotations")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
GenomeInfoDb::seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
class(annotation)
GenomeInfoDb::genome(annotation) <- "hg38"
class(annotation)

print("Setting up ATAC")
all_data[["ATAC"]] <- Signac::CreateChromatinAssay(
  counts = feature_mat,
  sep = c("-", "-"),
  fragments = fragments,
  annotation = annotation
)

save(date_of_run, session_info, 
     all_data,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2_part3.RData")

print("Basic processing of ATAC")
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data <- Signac::RunTFIDF(all_data)
all_data <- Signac::FindTopFeatures(all_data, min.cutoff = 'q0')
all_data <- Signac::RunSVD(all_data)

save(date_of_run, session_info, 
     all_data,
     file = "../../../../out/kevin/Writeup6l/Writeup6l_day0-macs2_part3.RData")

table(GenomicRanges::granges(all_data[["ATAC"]])@seqnames)

