rm(list=ls())
load("../../../../out/kevin/Writeup3d/09302021_sydney_basic_preprocess.RData")

gene_names <- rownames(all_data)
all_data <- Seurat::SCTransform(all_data,
                                residual.features = gene_names,
                                verbose = T)
dim(all_data[["SCT"]]@scale.data)

lin_mat <- all_data[["lineage"]]@counts
tmp <- Matrix::t(lin_mat)
lin_idx_list <- lapply(1:ncol(tmp), function(j){
  .nonzero_col(tmp, j)
})
names(lin_idx_list) <- colnames(tmp)
factor_vec <- as.factor(all_data@meta.data$Original_condition)
tabulate_mat <- t(sapply(lin_idx_list, function(idx){
  table(factor_vec[idx])
}))
rownames(tabulate_mat) <- colnames(tmp)
tabulate_mat <- tabulate_mat[order(tabulate_mat[,"naive"], decreasing = T),]
tmp <- tabulate_mat[,"cocl2"]/tabulate_mat[,"naive"]
tmp[is.infinite(tmp)] <- 0
lineage_ordering <- rownames(tabulate_mat)[order(tmp, decreasing = T)]
lineage_selected <-  rownames(tabulate_mat)[which(tmp > 10)]
length(lineage_selected)

# find all the relevant terminal cells in preparation for a differential expression
celltypes <- unique(all_data@meta.data$Original_condition)
cell_terminal <- which(all_data@meta.data$Original_condition == "cocl2")
tmp <- Matrix::t(lin_mat)
all_cell_names <- unique(unlist(lapply(lineage_selected, function(lineage){
  j <- which(colnames(tmp) == lineage)
  rownames(tmp)[.nonzero_col(tmp, j)]
})))
only_naive <- rownames(all_data@meta.data)[which(all_data@meta.data$Original_condition == "naive")]
naive_terminal <- all_cell_names[which(all_cell_names %in% only_naive)]
length(naive_terminal)

newgroup <- all_data@meta.data$Original_condition
names(newgroup) <- rownames(all_data@meta.data)
newgroup[naive_terminal] <- "naive_surviveCoCl2"
newgroup[setdiff(only_naive, naive_terminal)] <- "naive_noCoCl2"
all_data[["cocl2Status"]] <- newgroup
Seurat::Idents(all_data) <- "cocl2Status"
table(Seurat::Idents(all_data))

set.seed(10)
seurat_de <- Seurat::FindMarkers(all_data,
                                 assay = "SCT",
                                 slot = "scale.data",
                                 ident.1 = "naive_surviveCoCl2",
                                 ident.2 = "naive_noCoCl2",
                                 logfc.threshold = 0.1,
                                 min.pct = 0.05,
                                 verbose = T)
seurat_de[1:50,]
genes1 <- rownames(seurat_de)[which(seurat_de$p_val_adj <= 1e-3)]
genes2 <- rownames(seurat_de)[which(abs(seurat_de$avg_diff) > 0.1)]
genes <- intersect(genes1, genes2)
length(genes)

all_data[["SCT2"]] <- Seurat::CreateAssayObject(data = all_data[["SCT"]]@scale.data)
plot1 <- Seurat::DotPlot(all_data, features = genes[1:15],
                         assay = "SCT2",
                         group.by = "cocl2Status") + ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_sydney_basic_naivegenes.png"),
                plot1, device = "png", width = 12, height = 5, units = "in")

for(i in 1:15){
  plot1 <- Seurat::VlnPlot(all_data, features = genes[i],
                           idents = c("naive_surviveCoCl2", "naive_noCoCl2"),
                           assay = "SCT",
                           slot = "scale.data") + ggplot2::theme_classic()
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_sydney_basic_naivegenes_vln_", genes[i], ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

##############################

# custom list: 
genes_nancy <- c("DEAF1", "SERINC1", "NUDT8", "AP5B1", "CADM1", "PPP1R15A",
           "SERGEF", "NFATC2", "MARCH6", "GLS", "PRXL2A", "MT-ND6",
           "PREX1", "VPS53", "ROCK1")
seurat_de[genes_nancy,]

custom_mat <- all_data[["SCT"]]@scale.data[genes,only_naive]
custom_mat <- t(custom_mat)
set.seed(10)
custom_umap <- Seurat::RunUMAP(custom_mat)
custom_umap <- custom_umap@cell.embeddings
col_vec <- rep(rgb(0,0,0,0.2), nrow(custom_umap))
col_vec[which(rownames(custom_umap) %in% naive_terminal)] <- "red"
png("../../../../out/figures/Writeup3d/Writeup3d_sydney_basic_naivegenes_onlyumap.png",
    width = 1500, height = 1500, res = 300, units = "px")
plot(custom_umap[,1], custom_umap[,2], col = col_vec, pch = 16, asp = T,
     main = "UMAP, only Naive cells")
graphics.off()

#####

svd_res <- svd(custom_mat)
dimred <- svd_res$u[,1:2]%*%diag(svd_res$d[1:2])
col_vec <- rep(rgb(0,0,0,0.2), nrow(custom_umap))
col_vec[which(rownames(custom_umap) %in% naive_terminal)] <- "red"
png("../../../../out/figures/Writeup3d/Writeup3d_sydney_basic_naivegenes_onlypca.png",
    width = 1500, height = 1500, res = 300, units = "px")
plot(dimred[,1], dimred[,2], col = col_vec, pch = 16, asp = T,
     main = "PCA, only Naive cells")
graphics.off()


custom_mat <- all_data[["SCT"]]@scale.data[genes[1:5],only_naive]
custom_mat <- t(custom_mat)
col_vec <- rep(rgb(0,0,0,0.2), nrow(custom_umap))
col_vec[which(rownames(custom_umap) %in% naive_terminal)] <- "red"
png("../../../../out/figures/Writeup3d/Writeup3d_sydney_basic_naivegenes_pairs.png",
    width = 1500, height = 1500, res = 300, units = "px")
pairs(custom_mat, pch = 16, col = col_vec)
graphics.off()

