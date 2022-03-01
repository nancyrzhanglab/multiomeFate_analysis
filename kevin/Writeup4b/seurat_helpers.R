# from https://satijalab.org/signac/articles/pbmc_multiomic.html
# from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
# see https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/matrices for features.tsv annotation
create_seurat_object <- function(file_folder){
  tmp <- Seurat::Read10X_h5(paste0(file_prefix, file_folder, file_suffix))
  fragpath <- paste0(file_prefix, file_folder, "/outs/atac_fragments.tsv.gz")
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = tmp[["Gene Expression"]],
    assay = "RNA"
  )
  
  # create ATAC assay and add it to the object
  seurat_obj[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = tmp[["Peaks"]],
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = fragpath,
    annotation = annotation
  )
  
  seurat_obj
}

normalize_seurat <- function(seurat_obj){
  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 20000)
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj
}

qc_metrics <- function(seurat_obj){
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(object = seurat_obj, pattern = "^RPS")
  seurat_obj <- Seurat::CellCycleScoring(seurat_obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
  seurat_obj
}