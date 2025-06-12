rm(list=ls())
library(Seurat)
load("~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_subset.RData")

keep_vec <- rep(TRUE, ncol(seurat_object))
keep_vec[which(seurat_object$Cell.type.annotation == "Erythroid")] <- FALSE
seurat_object$keep <- keep_vec
seurat_object <- subset(seurat_object, keep == TRUE)

table(seurat_object$Time.point, seurat_object$Cell.type.annotation)

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
table(sapply(lineage_idx, length))

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

date_of_run <- Sys.time()
session_info <- devtools::session_info()
note <- paste("Putting in the assigned_lineage and keeping only cells with a lineage.")

save(seurat_object,
     date_of_run, session_info, note,
     file = "~/project/Multiome_fate/out/kevin/Writeup8/Writeup8_larry-dataset_step2_lineage-plotting.RData")

###############

length(unique(seurat_object$assigned_lineage))
tab_mat_full <- table(seurat_object$assigned_lineage, seurat_object$Time.point, seurat_object$Cell.type.annotation)
lineage_size <- apply(tab_mat_full, 1, sum)
tab_mat_full <- tab_mat_full[order(lineage_size, decreasing = T),,]

tab_mat <- table(seurat_object$assigned_lineage, seurat_object$Cell.type.annotation)
lineage_size <- apply(tab_mat, 1, sum)
tab_mat <- tab_mat[order(lineage_size, decreasing = T),]

df <- as.data.frame(apply(as.matrix.noquote(log10(tab_mat+1)),2,as.numeric))
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_step2_lineage-sizes_celltype.png",
                p1, device = "png", width = 8, height = 8, units = "in")

tab_mat <- table(seurat_object$assigned_lineage, seurat_object$Time.point)
lineage_size <- apply(tab_mat, 1, sum)
tab_mat <- tab_mat[order(lineage_size, decreasing = T),]

df <- as.data.frame(apply(as.matrix.noquote(log10(tab_mat+1)),2,as.numeric))
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_step2_lineage-sizes_timepoint.png",
                p1, device = "png", width = 8, height = 8, units = "in")

####################

tab_mat <- table(seurat_object$assigned_lineage, seurat_object$time_celltype)

df <- as.data.frame(apply(as.matrix.noquote(log10(tab_mat+1)),2,as.numeric))
p1 <- GGally::ggpairs(df, 
                      lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                      progress = FALSE) 
ggplot2::ggsave(filename = "~/project/Multiome_fate/out/figures/Writeup8/Writeup8_step2_lineage-sizes_time-celltype.png",
                p1, device = "png", width = 8, height = 8, units = "in")



