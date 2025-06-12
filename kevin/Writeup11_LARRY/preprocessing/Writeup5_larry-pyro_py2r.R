rm(list=ls()); gc(TRUE)
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(biomaRt)
options(Seurat.object.assay.version = "v5")

# from https://github.com/theislab/zellkonverter/issues/38
adata <- zellkonverter::readH5AD(
  file = "~/kzlinlab/data/larry_hematopoiesis_pyro-velocity/larry_kzlin-cleaned.h5ad"
)

names(adata@assays)

tmp <- Seurat::as.Seurat(adata, counts = "X", data = "X")
# converted most of the things. It creates an assay by default called "originalexp"
# see also https://www.bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html

seurat_obj <- Seurat::CreateSeuratObject(
  counts = tmp[["originalexp"]],
  data = tmp[["originalexp"]],
  meta.data = tmp@meta.data
)

seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")

# put in the assays
name_vec <- names(adata@assays)
gene_vec <- SeuratObject::Features(seurat_obj[["RNA"]])

# put in the gene metafeatures
gene_metadata <- SingleCellExperiment::rowData(adata)
seurat_obj[["RNA"]]@misc <- as.data.frame(gene_metadata)

# put in the dimension reductions
name_vec <- SingleCellExperiment::reducedDimNames(adata)
for(name_val in name_vec){
  mat <- SingleCellExperiment::reducedDim(adata, name_val)
  name_val2 <- paste0("python_", name_val)
  colnames(mat) <- paste0(name_val2, "_", 1:ncol(mat))

  seurat_obj[[name_val2]] <- Seurat::CreateDimReducObject(embeddings = mat,
                                                          assay = "RNA")
}

# put in the metadata
metadata_list <- adata@metadata
idx <- which(sapply(1:length(metadata_list), function(i){
  any(class(metadata_list[[i]]) %in% c("dgCMatrix", "dgRMatrix"))
}))
if(length(idx) > 0){
  graph_list <- metadata_list[idx]
  metadata_list <- metadata_list[-idx]
} else {
  graph_list <- numeric(0)
}
seurat_obj@misc <- metadata_list

if(length(graph_list) > 0){
  for(name_val in names(graph_list)){
    print(paste0("Putting in graph ", name_val))

    seurat_obj@graphs[[name_val]] <- graph_list[[name_val]]
    rownames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)
    colnames(seurat_obj@graphs[[name_val]]) <- SeuratObject::Cells(seurat_obj)

  }
}

Seurat::DefaultAssay(seurat_obj) <- "RNA"

# for the cluster color specifically, add the names for convenience
seurat_obj@misc[["state_info_colors"]] <- c(
  Baso = grDevices::rgb(245, 107, 0, maxColorValue = 255),
  Ccr7_DC = grDevices::rgb(6, 248, 162, maxColorValue = 255),
  Eos = grDevices::rgb(205, 1, 253, maxColorValue = 255),
  Erythroid = grDevices::rgb(71, 2, 250, maxColorValue = 255),
  Lymphoid = grDevices::rgb(0, 192, 255, maxColorValue = 255),
  Mast = grDevices::rgb(251, 0, 147, maxColorValue = 255),
  Meg = grDevices::rgb(6, 245, 38, maxColorValue = 255),
  Monocyte = grDevices::rgb(0, 55, 255, maxColorValue = 255),
  Neutrophil = grDevices::rgb(227, 172, 0, maxColorValue = 255),
  Undifferentiated = "gray"
)

# clearing environment
ls_vec <- ls()
ls_vec <- ls_vec[ls_vec != "seurat_obj"]
rm(list = ls_vec)
gc(TRUE)

################

# now do the usual seurat processing
Seurat::DefaultAssay(seurat_obj) <- "RNA"

# readjust the clonal modality
zz <- seurat_obj[["python_X_clone"]]@cell.embeddings
zz <- Matrix::Matrix(zz, sparse = TRUE)

## see https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
# if you want to find the nonzero entries for a row, I suggest
# first transposing via Matrix::t()
.nonzero_col <- function(mat, col_idx, bool_value){
  stopifnot(inherits(mat, "dgCMatrix"), col_idx %% 1 == 0,
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

lineage_list <- sapply(1:ncol(zz), function(j){
  .nonzero_col(mat = zz, col_idx = j, bool_value = FALSE)
})
stopifnot(as.numeric(names(table(table(unlist(lineage_list))))) == 1)
names(lineage_list) <- paste0("Lineage_", 1:length(lineage_list))
lineage_vec <- rep(NA, length(Seurat::Cells(seurat_obj)))
for(j in 1:length(lineage_list)){
  lineage_vec[lineage_list[[j]]] <- names(lineage_list[j])
}

seurat_obj$assigned_lineage <- lineage_vec
seurat_obj[["python_X_clone"]] <- NULL

#################

set.seed(10)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                           selection.method = "vst",
                                           nfeatures = 3000)
seurat_obj <- Seurat::ScaleData(seurat_obj)

set.seed(10)
seurat_obj <- Seurat::RunPCA(seurat_obj,
                             assay = "RNA",
                             features = Seurat::VariableFeatures(seurat_obj),
                             verbose = FALSE)
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              dims = 1:30)

####################

pdf(paste0("~/kzlinlab/projects/genePathwayManifold/git/genePathwayManifold_kevin/figs/kevin/Writeup5/Writeup5_larry-pyro_umaps.pdf"),
    onefile = T, width = 7, height = 5)

Seurat::DimPlot(seurat_obj,
                group.by = "state_info",
                cols = seurat_obj@misc[["state_info_colors"]],
                reduction = "python_X_emb")

Seurat::DimPlot(seurat_obj,
                group.by = "state_info",
                cols = seurat_obj@misc[["state_info_colors"]],
                reduction = "umap")

Seurat::DimPlot(seurat_obj,
                group.by = "time_info",
                reduction = "umap")

graphics.off()

SeuratObject::LayerData(seurat_obj,
                        assay = "RNA",
                        layer = "counts")[1:10,1:10]


mat <- SeuratObject::LayerData(seurat_obj,
                               assay = "RNA",
                               layer = "counts")
quantile(sapply(1:ncol(mat), function(j){length(.nonzero_col(mat,
                                                             col_idx = j,
                                                             bool_value = FALSE))/nrow(mat)}))

date_of_run <- Sys.time()
session_info <- devtools::session_info()

save(seurat_obj,
     date_of_run, session_info,
     file = "~/kzlinlab/projects/genePathwayManifold/out/kevin/Writeup5/Writeup5_larry-pyro.RData")

print("Done! :)")
