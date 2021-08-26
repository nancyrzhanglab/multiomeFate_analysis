rm(list=ls())

library(Seurat); library(Signac)

# load the original data with all the ATAC peak information
mbrain <- readRDS("../../../../data/Embrian_mouse/obj_seurat.rds")
# port over labels
mbrain_jane <- readRDS("../../../../data/Embrian_mouse/data_tenx_labels_Jane.rds")
cell_id <- rownames(mbrain@meta.data)
idx <- sapply(1:length(cell_id), function(x){
  which(rownames(mbrain_jane@meta.data) == cell_id[x])[1]
})
mbrain@meta.data$celltype <- mbrain_jane$savercatLable[idx]

load("../../../../out/kevin/Writeup3c/10x_mbrain_subset.RData")

##############
# select only the cells we want

cell_idx <- which(!mbrain2@meta.data$seurat_clusters %in% c(5, 7, 10, 12, 13, 14))
cell_id <- rownames(mbrain2@meta.data)[cell_idx]
keep_indicator <- as.numeric(sapply(rownames(mbrain@meta.data), function(x){x %in% cell_id}))
mbrain[["keep"]] <- keep_indicator
seurat_clusters <- rep(-1, nrow(mbrain@meta.data))
for(i in 1:nrow(mbrain2@meta.data)){
  cell_id <- rownames(mbrain2@meta.data)[i]
  idx <- which(rownames(mbrain@meta.data) == cell_id)
  if(length(idx) == 0) next()
  seurat_clusters[idx[1]] <- as.numeric(as.character(mbrain2@meta.data$seurat_clusters[i]))
}
mbrain[["new_seurat_clusters"]] <- seurat_clusters
mbrain3 <- subset(x = mbrain, keep == 1)

# now find the markers
n <- nrow(mbrain3@meta.data)
terminal_list <- list(16, 6, c(1,2,4), 9)
cell_terminal_list <- lapply(terminal_list, function(x){
  rownames(mbrain3@meta.data)[which(mbrain3@meta.data$new_seurat_clusters %in% x)]
})
initial_celltype <- 15
cell_initial_vec <- rownames(mbrain3@meta.data)[which(mbrain3@meta.data$new_seurat_clusters %in% initial_celltype)]

Seurat::DefaultAssay(mbrain3) <- "SCT"
de_list <- lapply(1:length(terminal_list), function(i){
  set.seed(10)
  Seurat:::FindMarkers(mbrain3[["SCT"]], cells.1 = cell_terminal_list[[i]],
                       cells.2 = unlist(cell_terminal_list[-i]), only.pos = T)
})
names(de_list) <- c("Oligo_vs_terminal", "Forebrain_vs_terminal",
                    "Cortical1_vs_terminal", "Cortical2_vs_terminal")
de_list2 <- lapply(1:length(terminal_list), function(i){
  set.seed(10)
  Seurat:::FindMarkers(mbrain3[["SCT"]], cells.1 = cell_terminal_list[[i]],
                       cells.2 = cell_initial_vec, only.pos = T)
})
names(de_list2) <- c("Oligo_vs_radial", "Forebrain_vs_radial",
                    "Cortical1_vs_radial", "Cortical2_vs_radial")
set.seed(10)
de_res3 <- Seurat:::FindMarkers(mbrain3[["SCT"]], cells.1 = cell_initial_vec,
                                 cells.2 = unlist(cell_terminal_list), only.pos = T)
de_combined <- c(de_list, de_list2)
de_combined[[length(de_combined)+1]] <- de_res3
names(de_combined)[length(de_combined)] <- "Radial_vs_terminal"

gene_vec <- unlist(lapply(de_combined, function(obj){
  idx <- which(obj$avg_log2FC > 0)
  rownames(obj[idx[1:min(length(idx), 50)],])
}))
table(table(gene_vec))
gene_vec <- sort(unique(gene_vec))

##############################

# now find the linked peaks to said genes
# [[note to self: is this the correct BSgenome to use?]]
library(BSgenome.Mmusculus.UCSC.mm10)
Seurat::DefaultAssay(mbrain3) <- "ATAC"
mbrain3 <- Signac::RegionStats(mbrain3, genome = BSgenome.Mmusculus.UCSC.mm10)

set.seed(10)
mbrain3 <- Signac::LinkPeaks(
  object = mbrain3,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = gene_vec
)

save(mbrain3, de_combined, file = "../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed_de.RData")


