rm(list=ls())
library(Seurat)

##### Create seurat obj #####

# download files from https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation
data_folder <- "~/nzhanglab/data/AllonKlein_hematopoietic_diff/"
out_folder <- "~/project/Multiome_fate/out/kevin/Writeup11b/"

# load in the metadata
cell_metadata_df <- read.csv(paste0(data_folder, "stateFate_inVitro_metadata.txt.gz"),
                             sep = "\t")
rowname_vec <- sapply(1:nrow(cell_metadata_df), function(i){
  paste0(cell_metadata_df[i,"Cell.barcode"], "_", i)
})
rownames(cell_metadata_df) <- rowname_vec

pseudotime_df <- read.csv(paste0(data_folder, "stateFate_inVitro_neutrophil_pseudotime.txt.gz"),
                          sep = "\t")
pseudotime_vec <- rep(NA, length(rowname_vec))
pseudotime_vec[1+pseudotime_df$Cell.index] <- pseudotime_df$pseudotime
trajectory_df <- read.csv(paste0(data_folder, "stateFate_inVitro_neutrophil_monocyte_trajectory.txt.gz"),
                          sep = "\t")
trajectory_vec <- rep(FALSE, length(rowname_vec))
trajectory_vec[trajectory_df$Cell.index+1] <- TRUE

cell_metadata_df$pseudotime <- pseudotime_vec
cell_metadata_df$trajectory <- trajectory_vec

# create the count matrix
expression_matrix <- Seurat::ReadMtx(mtx = paste0(data_folder, "stateFate_inVitro_normed_counts.mtx.gz"),
                                     cells = paste0(data_folder, "stateFate_inVitro_metadata.txt.gz"),
                                     features = paste0(data_folder, "stateFate_inVitro_gene_names.txt.gz"),
                                     cell.column = 2,
                                     feature.column = 1,
                                     mtx.transpose = TRUE,
                                     skip.cell = 1)
colnames(expression_matrix) <- rowname_vec
colsum_vec <- Matrix::colSums(expression_matrix)
stopifnot(diff(range(colsum_vec)) <= 1e-4)
expression_matrix <- expression_matrix * 1e6/colsum_vec[1]
expression_matrix@x <- round(expression_matrix@x)

# create the barcode matrix
barcode_mat <- Matrix::readMM(paste0(data_folder, "stateFate_inVitro_clone_matrix.mtx.gz"))
barcode_mat <- Matrix::t(barcode_mat)
colnames(barcode_mat) <- rowname_vec
rownames(barcode_mat) <- paste0("Lineage_", 1:nrow(barcode_mat))


#########################

seurat_obj <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                         data = expression_matrix,
                                         meta.data = cell_metadata_df)
seurat_obj[["Lineage"]] <- Seurat::CreateAssayObject(counts = barcode_mat)

Seurat::DefaultAssay(seurat_obj) <- "RNA"
set.seed(10)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                           selection.method = "vst",
                                           nfeatures = 2000)
seurat_obj <- Seurat::ScaleData(seurat_obj)
seurat_obj <- Seurat::RunPCA(seurat_obj,
                             features = Seurat::VariableFeatures(object = seurat_obj),
                             verbose = FALSE)
seurat_obj <- Seurat::RunUMAP(seurat_obj,
                              dims = 1:30)

# also include the author-generated SPRING-embedding
spring_embedding <- cbind(spring_1 = seurat_obj$SPRING.x,
                          spring_2 = seurat_obj$SPRING.y)
seurat_obj[["SPRING"]] <- Seurat::CreateDimReducObject(embeddings = spring_embedding)

# remove all the cells not in the trajectory
seurat_obj <- subset(seurat_obj, trajectory == TRUE)
seurat_obj <- subset(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))

mat <- SeuratObject::LayerData(object = seurat_obj, 
                               assay = "RNA", 
                               layer = "counts")
mat[201:210,1:10]

###########

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Basic formation of Seurat object from the LARRY dataset in https://github.com/AllonKleinLab/paper-data/tree/master/Lineage_tracing_on_transcriptional_landscapes_links_state_to_fate_during_differentiation")

save(seurat_obj,
     date_of_run, session_info, note,
     file = paste0(out_folder, "Writeup11b_larry-normalized.RData"))

print("Done! :)")
