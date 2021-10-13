rm(list=ls())
library(Seurat)
load("../../../../out/kevin/Writeup3d/09302021_sydney_basic_preprocess.RData")

gene_names <- rownames(all_data)
all_data <- Seurat::SCTransform(all_data,
                                residual.features = gene_names,
                                do.center = F,
                                verbose = T)
dim(all_data[["SCT"]]@scale.data)

save(all_data, 
     file = "../../../../out/kevin/Writeup3e/10122021_sydney_sctransform_preprocess.RData")

###############

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
head(tabulate_mat)
quantile(tabulate_mat[,"naive"], probs = seq(0.9,1,length.out=11))
tmp <- tabulate_mat[,"cocl2"]/tabulate_mat[,"naive"]
tmp[is.infinite(tmp)] <- 0
lineage_idx <- intersect(which(tmp > 1), which(tabulate_mat[,"naive"] >= 50))
length(lineage_idx)
lineage_selected <-  rownames(tabulate_mat)[lineage_idx]
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
length(only_naive)

##############################

# set up identities
newgroup <- all_data@meta.data$Original_condition
names(newgroup) <- rownames(all_data@meta.data)
newgroup[naive_terminal] <- "naive_surviveCoCl2"
newgroup[setdiff(only_naive, naive_terminal)] <- "naive_noCoCl2"
all_data[["cocl2Status"]] <- newgroup
Seurat::Idents(all_data) <- "Original_condition"
table(Seurat::Idents(all_data))

##############################

# first make sure we can separate cocl2 from naive
marker_genes <- c("AXL", "EGFR", "JUN", "VEGFC", "WNT5A", "NGFR",
                  "SERPINE1", "FGFR1", "LOXL2", "NRG1", "PDGFRB")
set.seed(10)
seurat_de <- Seurat::FindMarkers(all_data,
                                 assay = "SCT",
                                 slot = "scale.data",
                                 ident.1 = "cocl2",
                                 ident.2 = "naive",
                                 logfc.threshold = 0.1,
                                 min.pct = 0.05,
                                 verbose = T)
seurat_de[c("AXL", "NRG1", "EGFR"),]
seurat_de[1:50,]
head(seurat_de)

for(i in 1:length(marker_genes)){
  plot1 <- Seurat::VlnPlot(all_data, features = marker_genes[i],
                           idents = c("cocl2", "naive"),
                           assay = "SCT",
                           slot = "scale.data") + ggplot2::theme_classic()
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/Writeup3e_sydney_Cocl2vsNaive_vln_", marker_genes[i], ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

mean(all_data[["SCT"]]@scale.data["EGFR",])
############################

# now do the among-naive comparisons
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

tentative_genes <- rownames(seurat_de)[which(seurat_de$p_val_adj <= 1e-3)]
length(tentative_genes)
idx <- which(all_data@meta.data$Original_condition == "naive")
mat <- all_data[["SCT"]]@scale.data[tentative_genes,idx]

# see if there's a clear separation between low and high thresholds
cells.1 <- which(colnames(mat) %in% naive_terminal)
cells.2 <- setdiff(1:ncol(mat), cells.1)
proportion_mat <- lapply(1:nrow(mat), function(j){
  vec <- mat[j,]
  min_val <- max(vec[vec <= 0])
  max_val <- min(vec[vec >= 0])
  sd1 <- stats::sd(vec[vec <= 0])
  sd2 <- stats::sd(vec[vec > 0])
  sd_val <- min(c(sd1, sd2))
  
  bool <- abs(min_val - max_val) >= 3*sd_val
  if(!bool || sd1 > sd2) return(list(mat = NA, pval = 1))
  
  res_mat <- matrix(0, 2, 2)
  colnames(res_mat) <- c("survive", "die")
  rownames(res_mat) <- c("low", "high")
  res_mat[1,1] <- length(which(vec[cells.1] <= 0))
  res_mat[2,1] <- length(which(vec[cells.1] > 0))
  
  res_mat[1,2] <- length(which(vec[cells.2] <= 0))
  res_mat[2,2] <- length(which(vec[cells.2] > 0))
  
  prop_test <- stats::prop.test(x = res_mat[2,], n = colSums(res_mat))
  
  list(mat = res_mat, pval = prop_test$p.value)
})
names(proportion_mat) <- rownames(mat)
marker_genes <- names(proportion_mat)[order(sapply(proportion_mat, function(x){
  if(!all(is.na(x$mat))){
    max(x$mat[2,1]/x$mat[2,2], x$mat[2,2]/x$mat[2,1])
  } else {
    0
  }
}), decreasing = T)[1:15]]

marker_genes <- unique(c(marker_genes, "AXL", "EGFR", "NRG1"))
for(i in 1:length(marker_genes)){
  plot1 <- Seurat::VlnPlot(all_data, features = marker_genes[i],
                           idents = c("naive_surviveCoCl2", "naive_noCoCl2"),
                           assay = "SCT",
                           slot = "scale.data") + ggplot2::theme_classic()
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup3e/Writeup3e_sydney_naive_withinCocl2_vln_", marker_genes[i], ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}
