rm(list=ls())
library(Seurat)
##### Create seurat obj #####

# download files from https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation

# load in the metadata
cell_metadata_df <- read.csv("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_metadata.txt.gz",
                             sep = "\t")
rowname_vec <- sapply(1:nrow(cell_metadata_df), function(i){
  paste0(cell_metadata_df[i,"Cell.barcode"], "_", i)
})
rownames(cell_metadata_df) <- rowname_vec

pseudotime_df <- read.csv("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_neutrophil_pseudotime.txt.gz",
                          sep = "\t")
pseudotime_vec <- rep(NA, length(rowname_vec))
pseudotime_vec[1+pseudotime_df$Cell.index] <- pseudotime_df$pseudotime
trajectory_df <- read.csv("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_neutrophil_monocyte_trajectory.txt.gz",
                          sep = "\t")
trajectory_vec <- rep(FALSE, length(rowname_vec))
trajectory_vec[trajectory_df$Cell.index+1] <- TRUE

cell_metadata_df$pseudotime <- pseudotime_vec
cell_metadata_df$trajectory <- trajectory_vec

# create the count matrix
expression_matrix <- Seurat::ReadMtx(mtx = "~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_normed_counts.mtx.gz",
                                     cells = "~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_metadata.txt.gz",
                                     features = "~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_gene_names.txt.gz",
                                     cell.column = 2,
                                     feature.column = 1,
                                     mtx.transpose = T,
                                     skip.cell = 1)
colnames(expression_matrix) <- rowname_vec

# create the barcode matrix
barcode_mat <- Matrix::readMM("~/nzhanglab/data/AllonKlein_hematopoietic_diff/stateFate_inVitro_clone_matrix.mtx.gz")
barcode_mat <- Matrix::t(barcode_mat)
colnames(barcode_mat) <- rowname_vec
rownames(barcode_mat) <- paste0("Lineage_", 1:nrow(barcode_mat))

#########################

seurat_object <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                            data = expression_matrix,
                                            meta.data = cell_metadata_df)
seurat_object[["Lineage"]] <- Seurat::CreateAssayObject(counts = barcode_mat)

Seurat::DefaultAssay(seurat_object) <- "RNA"
set.seed(10)
seurat_object <- Seurat::FindVariableFeatures(seurat_object,
                                              selection.method = "vst",
                                              nfeatures = 2000)

# also include the author-generated SPRING-embedding
spring_embedding <- cbind(spring_1 = seurat_object$SPRING.x,
                          spring_2 = seurat_object$SPRING.y)
seurat_object[["SPRING"]] <- Seurat::CreateDimReducObject(embeddings = spring_embedding)

###########

# remove all the cells not in the trajectory
seurat_object <- subset(seurat_object, trajectory == TRUE)

###########

keep_vec <- rep(TRUE, ncol(seurat_object))
keep_vec[which(seurat_object$Cell.type.annotation == "Erythroid")] <- FALSE
seurat_object$keep <- keep_vec
seurat_object <- subset(seurat_object, keep == TRUE)

# assigning cells to a lineage
mat <- SeuratObject::LayerData(seurat_object, layer = "counts", assay = "Lineage")
total_count <- Matrix::colSums(mat)
table(total_count)
n <- ncol(mat)
lineage_idx <- lapply(1:n, function(i){
  multiomeFate:::.nonzero_col(mat = mat, 
                              col_idx = i, 
                              bool_value = F)
})

lineage_vec <- rep(NA, ncol(seurat_object))
names(lineage_vec) <- SeuratObject::Cells(seurat_object)
for(i in 1:n){
  if(length(lineage_idx[[i]]) == 1){
    lineage_vec[i] <- rownames(mat)[lineage_idx[[i]]]
  }
}
seurat_object$assigned_lineage <- lineage_vec

# remove all the cells without a lineage
keep_vec <- rep(TRUE, ncol(seurat_object))
keep_vec[which(is.na(seurat_object$assigned_lineage))] <- FALSE
seurat_object$keep <- keep_vec
seurat_object <- subset(seurat_object, keep == TRUE)

time_celltype_df <- seurat_object@meta.data[,c("Cell.type.annotation", "Time.point")]
time_celltype_str <- apply(time_celltype_df, 1, function(x){paste0(x, collapse = "-")})
time_celltype_str <- factor(time_celltype_str)
seurat_object$time_celltype <- time_celltype_str

###########

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Basic formation of Seurat object from the LARRY dataset in https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation")

save(seurat_object,
     date_of_run, session_info, note,
     file = "~/project/Multiome_fate/out/kevin/Writeup13/Writeup13_larry-dataset.RData")

print("Done! :)")

