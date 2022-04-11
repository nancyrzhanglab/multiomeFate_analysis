rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

create_seurat_object <- function(file_prefix, file_folder, file_suffix){
  tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  fragpath <- paste0(file_prefix, file_folder, "/outs/atac_fragments.tsv.gz")
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = tmp[["Gene Expression"]],
    assay = "RNA"
  )
  
  seurat_obj
}

##########################

file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_04_arc_time0", "2022_04_arc_time10_CIS", 
                  "2022_04_arc_time10_COCL2", "2022_04_arc_time10_DABTRAM",
                  "2022_04_arc_week5_CIS", "2022_04_arc_week5_COCL2",
                  "2022_04_arc_week5_DABTRAM",
                  "2022_04_arc_test2", "2022_04_arc_test3")
name_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM",
              "week5_CIS", "week5_COCL2", "week5_DABTRAM", "test2", "test3")

annotation <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotation) <- "UCSC"
Signac::genome(annotation) <- "hg38"

seurat_list <- lapply(1:length(file_folders), function(i){
  file_folder <- file_folders[i]
  seurat_obj <- create_seurat_object(file_prefix = file_prefix,
                                     file_folder = file_folder,
                                     file_suffix = file_suffix)
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})

all_data <- merge(seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]], 
                                          seurat_list[[4]], seurat_list[[5]], 
                                          seurat_list[[6]], seurat_list[[7]],
                                          seurat_list[[8]], seurat_list[[9]]), 
                  add.cell.ids = name_vec, 
                  project = "All_Data", merge.data = T)

print("Finished processing RNA")
save(all_data, seurat_list, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4c/Writeup4c_timeAll_peakmerging.RData")

#######################

# https://satijalab.org/signac/articles/merging.html

# read in peak sets
peaks_list <- lapply(file_folders, function(file_folder){
  read.table(
    file = paste0(file_prefix, file_folder, "/outs/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
})

# convert to genomic ranges
gr_list <- lapply(peaks_list, function(peaks_data){
  GenomicRanges::makeGRangesFromDataFrame(peaks_data)
})

# Create a unified set of peaks to quantify in each dataset
combined_peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], 
                               gr_list[[3]], gr_list[[4]], 
                               gr_list[[5]], gr_list[[6]], 
                               gr_list[[7]], gr_list[[8]], 
                               gr_list[[9]]))

# Filter out bad peaks based on length
peak_widths <- IRanges::width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths < 10000 & peak_widths > 20]

# create fragment objects
frags_list <- lapply(1:length(file_folders), function(i){
  file_folder <- file_folders[i]
  Signac::CreateFragmentObject(
    path = paste0(file_prefix, file_folder, "/outs/atac_fragments.tsv.gz"),
    cells = colnames(seurat_list[[i]])
  )
})

atac_count_list <- lapply(1:length(file_folders), function(i){
  Signac::FeatureMatrix(
    fragments = frags_list[[i]],
    features = combined_peaks,
    cells = colnames(seurat_list[[i]])
  )
})

seurat_atac_list <- lapply(1:length(file_folders), function(i){
  assay_obj <- Signac::CreateChromatinAssay(counts = atac_count_list[[i]],
                                            fragments = frags_list[[i]],
                                            sep = c("-", "-"),
                                            annotation = annotation)
  seurat_obj <- Seurat::CreateSeuratObject(assay_obj, assay = "ATAC")
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})

all_data_atac <- merge(seurat_atac_list[[1]], y = c(seurat_atac_list[[2]], seurat_atac_list[[3]], 
                                                    seurat_atac_list[[4]], seurat_atac_list[[5]], 
                                                    seurat_atac_list[[6]], seurat_atac_list[[7]],
                                                    seurat_atac_list[[8]], seurat_atac_list[[9]]), 
                  add.cell.ids = name_vec, 
                  project = "All_Data", merge.data = T)

print("Finished processing ATAC")
save(all_data, all_data_atac, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4c/Writeup4c_timeAll_peakmerging.RData")

###################

all_data[["ATAC"]] <- all_data_atac[["ATAC"]]

Seurat::DefaultAssay(all_data) <- "ATAC"
all_data <- Signac::NucleosomeSignal(object = all_data)
all_data <- Signac::TSSEnrichment(all_data, fast = FALSE)
all_data$high.tss <- ifelse(all_data$TSS.enrichment > 2, 'High', 'Low')
all_data$blacklist_fraction <- Signac::FractionCountsInRegion(
  object = all_data, 
  assay = "ATAC",
  regions = blacklist_hg19
)

meta_list <- lapply(file_folders, function(file_folder){
  read.table(
    file = paste0(file_prefix, file_folder, "/outs/per_barcode_metrics.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
})

all_meta <- Reduce(rbind, meta_list)
passed_filter_vec <- rep(NA, ncol(all_data))
names(passed_filter_vec) <- colnames(all_data)
peak_region_fragments_vec <- rep(NA, ncol(all_data))
names(peak_region_fragments_vec) <- colnames(all_data)
tmp_name <- sapply(names(passed_filter_vec), function(name_val){
  tmp <- strsplit(name_val, split = "_")[[1]]
  tmp[length(tmp)]
})
names(tmp_name) <- NULL
all_meta <- all_meta[which(rownames(all_meta) %in% tmp_name), ]
for(i in 1:length(passed_filter_vec)){
  if(i %% floor(length(passed_filter_vec)/10) == 0) cat('*')
  passed_filter_vec[i] <- all_meta[tmp_name[i],"atac_fragments"]
  peak_region_fragments_vec[i] <- all_meta[tmp_name[i],"atac_peak_region_fragments"]
}
all_data$passed_filters <- passed_filter_vec
all_data$peak_region_fragments <- peak_region_fragments_vec
all_data$pct_reads_in_peaks <- all_data$peak_region_fragments / all_data$passed_filters * 100

Seurat::DefaultAssay(all_data) <- "RNA"
all_data[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = all_data, pattern = "^MT-")
all_data[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = all_data, pattern = "^RPS")

print("Finished processing metadata")
save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4c/Writeup4c_timeAll_peakmerging.RData")

##########################

Seurat::DefaultAssay(all_data) <- "RNA"
set.seed(10)
all_data <- Seurat::FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 500)
jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC")
jackpot_genes <- intersect(jackpot_genes, rownames(all_data))
all_data[["RNA"]]@var.features <- unique(c(jackpot_genes, all_data[["RNA"]]@var.features))
all_data <- Seurat::ScaleData(all_data)
all_data <- Seurat::RunPCA(all_data, verbose = FALSE) 
set.seed(10)
all_data <- Seurat::RunUMAP(all_data, dims = 1:50)
all_data <- Seurat::CellCycleScoring(all_data, 
                                     g2m.features = cc.genes$g2m.genes, 
                                     s.features = cc.genes$s.genes)

Seurat::DefaultAssay(all_data) <- "ATAC"
set.seed(10)
all_data <- Signac::RunTFIDF(all_data)
all_data <- Signac::FindTopFeatures(all_data, min.cutoff = 'q0')
all_data <- Signac::RunSVD(all_data)
set.seed(10)
all_data <- Seurat::RunUMAP(object = all_data, 
                            reduction = 'lsi', dims = 2:50,
                            reduction.name = 'adt.umap')

##################

dataset_vec <- sapply(rownames(all_data@meta.data), function(x){
  tmp <- strsplit(x, split = "_")[[1]]
  if(length(tmp) == 2){
    return(tmp[1])
  } else {
    return(paste0(tmp[1:2], collapse = "_"))
  }
})
names(dataset_vec) <- NULL
all_data$original_dataset <- dataset_vec

print("Finished preprocessing data")
save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4c/Writeup4c_timeAll_peakmerging.RData")

