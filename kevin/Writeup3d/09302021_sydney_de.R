rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/kevin/Writeup3d/09302021_sydney_de_preprocess.RData")
ls()

# tabulate the lineages
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

# find all the relevant terminal cells in preparation for a differential expression
celltypes <- unique(all_data@meta.data$Original_condition)
cell_terminal <- which(all_data@meta.data$Original_condition == "cocl2")
lineage_terminal <- which(sparseMatrixStats::rowSums2(lin_mat[,cell_terminal,drop=F]) > 0)
tmp <- Matrix::t(lin_mat)
all_cell_names <- unique(unlist(lapply(lineage_terminal, function(j){
  colnames(lin_mat)[.nonzero_col(tmp, j)]
})))
only_naive <- rownames(all_data@meta.data)[which(all_data@meta.data$Original_condition == "naive")]
naive_terminal <- all_cell_names[which(all_cell_names %in% only_naive)]
length(naive_terminal)

# cells.1 <- naive_terminal
# cells.2 <- setdiff(only_naive, naive_terminal)

newgroup <- all_data@meta.data$Original_condition
names(newgroup) <- rownames(all_data@meta.data)
newgroup[naive_terminal] <- "naive_surviveCoCl2"
newgroup[setdiff(only_naive, naive_terminal)] <- "naive_noCoCl2"
all_data[["cocl2Status"]] <- newgroup
Seurat::Idents(all_data) <- "cocl2Status"

set.seed(10)
seurat_de <- Seurat::FindMarkers(all_data,
                                 ident.1 = "naive_surviveCoCl2",
                                 ident.2 = "naive_noCoCl2",
                                 test.use = "negbinom",
                                 logfc.threshold = 0.1,
                                 min.pct = 0.05,
                                 verbose = T)
genes <- rownames(seurat_de)[which(seurat_de$p_val_adj <= 1e-2)]

plot1 <- Seurat::DotPlot(all_data, features = genes,
                         group.by = "cocl2Status") + ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_sydney_de_naivegenes.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

for(i in 1:length(genes)){
  plot1 <- Seurat::VlnPlot(all_data, features = genes[i],
                           idents = c("naive_surviveCoCl2", "naive_noCoCl2"),
                           assay = "SCT",
                           slot = "data") + ggplot2::theme_classic()
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3d/Writeup3d_sydney_de_naivegenes_vln_", genes[i], ".png"),
                  plot1, device = "png", width = 10, height = 5, units = "in")
}

