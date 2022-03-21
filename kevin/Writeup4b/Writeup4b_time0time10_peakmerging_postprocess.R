rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_peakmerging.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

all_data[["ATAC"]] <- all_data_atac[["ATAC"]]

Seurat::DefaultAssay(all_data) <- "ATAC"
all_data <- Signac::NucleosomeSignal(object = all_data)
plot1 <- Signac::FragmentHistogram(object = all_data, 
                                   region = 'chr1-1-10000000')
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_fragmentSignal_nucleosomeSignal.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

all_data <- Signac::TSSEnrichment(all_data, fast = FALSE)
all_data$high.tss <- ifelse(all_data$TSS.enrichment > 2, 'High', 'Low')
plot1 <- Signac::TSSPlot(all_data, group.by = 'high.tss') + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_tssEnrichment.png"),
                plot1, device = "png", width = 8, height = 4, units = "in")

# from https://satijalab.org/signac/articles/pbmc_vignette.html
all_data$blacklist_fraction <- Signac::FractionCountsInRegion(
  object = all_data, 
  assay = "ATAC",
  regions = blacklist_hg19
)

# see https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics
time0_meta <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time0/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
time10_cis_meta <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_CIS/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
time10_cocl2_meta <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_COCL2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
time10_dabtram_meta <- read.table(
  file = "../../../../BarcodeOutputs/2022_02/Cellranger_count_output/2022_02_arc_time10_DABTRAM/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
all_meta <- rbind(time0_meta, time10_cis_meta, time10_cocl2_meta, time10_dabtram_meta)

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

plot1 <- Seurat::VlnPlot(
  object = all_data,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
  group.by = "dataset", 
  pt.size = 0.1,
  ncol = 5
)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_atac_QC.png"),
                plot1, device = "png", width = 15, height = 4, units = "in")


save(all_data, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_processed.RData")

########################

rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
load("../../../../out/kevin/Writeup4b/Writeup4b_time0time10_alldata_GEXandATAC_processed.RData")

set.seed(10)
all_data <- Signac::RunTFIDF(all_data)
all_data <- Signac::FindTopFeatures(all_data, min.cutoff = 'q0')
all_data <- Signac::RunSVD(all_data)

set.seed(10)
all_data <- Seurat::RunUMAP(object = all_data, 
                            reduction = 'lsi', dims = 2:30)

plot1 <- Seurat::DimPlot(all_data, 
                         reduction = "umap", 
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (TF-IDF),\nusing 29 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_atac_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#####################

set.seed(10)
all_data <- Seurat::RunUMAP(object = all_data, 
                            reduction = 'lsi', dims = 2:50)
plot1 <- Seurat::DimPlot(all_data, 
                         reduction = "umap", 
                         group.by = "dataset", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("ATAC (TF-IDF),\nusing 49 PCs"))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4b/Writeup4b_atac_umap_50pc.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
