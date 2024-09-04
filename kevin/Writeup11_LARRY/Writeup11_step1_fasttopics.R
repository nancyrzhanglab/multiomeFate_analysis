rm(list=ls())
library(Seurat)
library(fastTopics)
load("~/nzhanglab/data/LARRY_pyro-velocity/Writeup5_larry-pyro.RData")

out_folder <- "~/project/Multiome_fate/out/kevin/Writeup11/"

any(duplicated(seurat_obj$Cell.barcode))

# put the trajectory information in
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

# line up the two datasets
tmp <- rownames(cell_metadata_df)
tmp <- sapply(tmp, function(x){strsplit(x, split = "_")[[1]][1]})
tmp <- gsub(pattern = "-", replace = "", tmp)
for(i in 1:length(tmp)){
  tmp[i] <- paste0(cell_metadata_df[i,"Library"], ":", tmp[i])
}
length(intersect(tmp, Seurat::Cells(seurat_obj)))
any(duplicated(tmp))
rownames(cell_metadata_df) <- tmp
cell_metadata_df <- cell_metadata_df[Seurat::Cells(seurat_obj),]
seurat_obj$pseudotime <- cell_metadata_df$pseudotime
seurat_obj$trajectory <- cell_metadata_df$trajectory

# remove all the cells not in the trajectory
seurat_obj <- subset(seurat_obj, trajectory == TRUE)

save(seurat_obj, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11_larry_seurat.RData"))

#########################

mat <- SeuratObject::LayerData(object = seurat_obj, 
                               assay = "RNA", 
                               layer = "counts",
                               features = Seurat::VariableFeatures(seurat_obj, assay = "RNA"))
mat <- Matrix::t(mat)

K <- 30
set.seed(10)
topic_res <- fastTopics::fit_topic_model(mat, k = K)

save(topic_res, date_of_run, session_info,
     file = paste0(out_folder, "Writeup11_larry_fasttopics.RData"))

print("Done! :)")
