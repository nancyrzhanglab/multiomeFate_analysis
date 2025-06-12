rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

create_seurat_object <- function(file_prefix, file_folder, file_suffix){
  tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = tmp[["Gene Expression"]],
    assay = "RNA"
  )
  
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj
}

##########################

file_prefix <- "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/"
file_suffix <- "/outs/filtered_feature_bc_matrix.h5"
file_folders <- c("2022_05_19_arc_time0", "2022_05_19_arc_time10_CIS", 
                  "2022_05_19_arc_time10_COCL2", "2022_05_19_arc_time10_DABTRAM",
                  "2022_05_19_arc_week5_CIS", "2022_05_19_arc_week5_COCL2",
                  "2022_05_19_arc_week5_DABTRAM")
name_vec <- c("day0", "day10_CIS", "day10_COCL2", "day10_DABTRAM",
              "week5_CIS", "week5_COCL2", "week5_DABTRAM")

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
                                          seurat_list[[6]], seurat_list[[7]]), 
                  add.cell.ids = name_vec, 
                  project = "All_Data", merge.data = T)

print("Finished processing RNA")
save(all_data, seurat_list, date_of_run, session_info,
     file = "../../../../out/kevin/analysis_pipeline/step1_peakmerging.RData")

#######################

# https://satijalab.org/signac/articles/merging.html

# read in peak sets
print("Reading peak set")
peaks_list <- lapply(file_folders, function(file_folder){
  read.table(
    file = paste0(file_prefix, file_folder, "/outs/atac_peaks.bed"),
    col.names = c("chr", "start", "end")
  )
})

print("Creating genomic ranges")
# convert to genomic ranges
gr_list <- lapply(peaks_list, function(peaks_data){
  GenomicRanges::makeGRangesFromDataFrame(peaks_data)
})

print("Creating unified peak set")
# Create a unified set of peaks to quantify in each dataset
combined_peaks <- reduce(x = c(gr_list[[1]], gr_list[[2]], 
                               gr_list[[3]], gr_list[[4]], 
                               gr_list[[5]], gr_list[[6]], 
                               gr_list[[7]]))

# Filter out bad peaks based on length
peak_widths <- IRanges::width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths < 10000 & peak_widths > 20]

print("Creating frag objects")
# create fragment objects
frags_list <- lapply(1:length(file_folders), function(i){
  file_folder <- file_folders[i]
  Signac::CreateFragmentObject(
    path = paste0(file_prefix, file_folder, "/outs/atac_fragments.tsv.gz"),
    cells = colnames(seurat_list[[i]])
  )
})

print("Constructing count matrices")
atac_count_list <- lapply(1:length(file_folders), function(i){
  Signac::FeatureMatrix(
    fragments = frags_list[[i]],
    features = combined_peaks,
    cells = colnames(seurat_list[[i]])
  )
})

print("Constructing Seurat objects")
seurat_atac_list <- lapply(1:length(file_folders), function(i){
  assay_obj <- Signac::CreateChromatinAssay(counts = atac_count_list[[i]],
                                            fragments = frags_list[[i]],
                                            sep = c("-", "-"),
                                            annotation = annotation)
  seurat_obj <- Seurat::CreateSeuratObject(assay_obj, assay = "ATAC")
  seurat_obj$dataset <- name_vec[i]
  seurat_obj
})

print("Merging ATAC peaks")
all_data_atac <- merge(seurat_atac_list[[1]], y = c(seurat_atac_list[[2]], seurat_atac_list[[3]], 
                                                    seurat_atac_list[[4]], seurat_atac_list[[5]], 
                                                    seurat_atac_list[[6]], seurat_atac_list[[7]]), 
                       add.cell.ids = name_vec, 
                       project = "All_Data", merge.data = T)

print("Finished processing ATAC")
save(all_data, all_data_atac, date_of_run, session_info,
     file = "../../../../out/kevin/analysis_pipeline/step1_peakmerging.RData")

###################

print("Preprocessing ATAC")
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
     file = "../../../../out/kevin/analysis_pipeline/step1_peakmerging.RData")
