rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)

# see also https://github.com/nancyrzhanglab/multiomeFate_analysis/blob/emilia/emilia/task_identify_features_predictive_growth_and_lineage_specific/Run_chromVAR_Vierstra.R

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup10a/"
load(paste0(out_folder, "Writeup10a_ppStep6_qc-step2.RData"))

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

## see https://stuartlab.org/signac/articles/motif_vignette.html
## https://stuartlab.org/signac/articles/data_structures.html#the-motif-class
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# https://github.com/stuart-lab/signac/issues/486
main.chroms <- GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(GenomeInfoDb::seqnames(GenomicRanges::granges(all_data[["ATAC"]]))) %in% main.chroms)
all_data[["ATAC"]] <- subset(all_data[["ATAC"]], 
                             features = SeuratObject::Features(all_data[["ATAC"]])[keep.peaks])

Seurat::DefaultAssay(all_data) <- "ATAC"
all_data[["RNA"]] <- NULL
all_data[["Lineage"]] <- NULL
all_data[["pca"]] <- NULL
all_data[["umap"]] <- NULL

# Run the following for each treatment-timepoint pair
dataset_levels <- unique(all_data$dataset)

for(dataset in dataset_levels){
  print(paste0("Working on ", dataset))
  
  keep_vec <- rep(FALSE, length(Seurat::Cells(all_data)))
  keep_vec[which(all_data$dataset == dataset)] <- TRUE
  all_data$keep <- keep_vec
  all_data_subset <- subset(all_data, keep == TRUE)
  
  # the above line didnt' work for some reason, so we'll run it piece by piece below
  ## see https://github.com/stuart-lab/signac/blob/master/R/motifs.R
  print("Constructing motif matrix")
  motif.matrix <- Signac::CreateMotifMatrix(
    features = GenomicRanges::granges(all_data_subset[["ATAC"]]),
    pwm = pfm,
    genome = "hg38",
    use.counts = FALSE
  )
  
  print("Matching motifs")
  motif.positions <- motifmatchr::matchMotifs(
    pwms = pfm,
    subject = GenomicRanges::granges(all_data_subset[["ATAC"]]),
    out = "positions",
    genome = "hg38"
  )
  
  motif <- Signac::CreateMotifObject(
    data = motif.matrix,
    positions = motif.positions,
    pwm = pfm
  )
  
  all_data_subset[["ATAC"]] <- Seurat::SetAssayData(
    object = all_data_subset[["ATAC"]],
    slot = 'motifs',
    new.data = motif
  )
  
  print("Running ChromVar")
  all_data_subset <- Signac::RunChromVAR(
    object = all_data_subset,
    genome = "hg38"
  )
  
  Seurat::DefaultAssay(all_data_subset) <- "ATAC"
  empty_matrix <- Matrix::Matrix(0, 
                                 nrow = length(SeuratObject::Features(all_data_subset)),
                                 ncol = length(Seurat::Cells(all_data_subset)),
                                 sparse = TRUE)
  rownames(empty_matrix) <- SeuratObject::Features(all_data_subset)
  colnames(empty_matrix) <- Seurat::Cells(all_data_subset)
  SeuratObject::LayerData(all_data_subset,
                          assay = "ATAC",
                          layer = "counts") <- empty_matrix
  SeuratObject::LayerData(all_data_subset,
                          assay = "ATAC",
                          layer = "data") <- empty_matrix
  
  save(date_of_run, session_info, 
       all_data_subset,
       file = paste0(out_folder, "Writeup10a_ppStep7d_chromvar_", dataset, ".RData"))
}

print("Done! :)")
