rm(list=ls())

load("../../../../out/kevin/Writeup3c/10x_mbrain_fate_preprocessed_de.RData")
chiyun <- readRDS("../../../../out/kevin/Writeup3c/chiyun_08282021_oligo_linkpeaks.rds")
library(Seurat); library(Signac); library(multiomeFate)

set.seed(10)
date_of_run <- Sys.time()

mat_x <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))

genes <- unique(c(rownames(de_combined[[1]][1:50,]), rownames(de_combined[[5]][1:50,])))
genes <- genes[which(genes %in% colnames(mat_y))]

precursor_score <- lapply(genes, function(gene){
  idx <- which(mbrain3[["ATAC"]]@links$gene == gene)
  peak_names <- mbrain3[["ATAC"]]@links$peak[idx]
  vec1 <- rowSums(mat_x[,which(colnames(mat_x) %in% peak_names), drop = F])
  vec2 <- mat_y[,which(colnames(mat_y) == gene)]
  
  tmp <- cbind(vec1, vec2)
  colnames(tmp) <- c("precursor", "rna")
  tmp
})

precursor_mat <- do.call(cbind, precursor_score)
colnames(precursor_mat) <- paste0("precursor", 1:ncol(precursor_mat))
mbrain3[["precursor"]] <- Seurat::CreateDimReducObject(embedding = precursor_mat, 
                                                    key = "precursor_", assay = "RNA")

for(i in 1:length(genes)){
  print(i)
  
  gene <- genes[i]
  
  plot1 <- Seurat::FeaturePlot(mbrain3, features = paste0("precursor_", 2*(i-1)+2), reduction = "wnn.umap")
  plot1 <- plot1 + ggplot2::labs(title = paste0(gene, ": RNA"), x = "wnnUMAP_1", y = "wnnUMAP_2")

  plot2 <- Seurat::FeaturePlot(mbrain3, features = paste0("precursor_", 2*(i-1)+1), reduction = "wnn.umap")
  plot2 <- plot2 + ggplot2::labs(title = paste0(gene, ": Linked ATAC (Sum)"), x = "wnnUMAP_1", y = "wnnUMAP_2")
  
  plot3 <- Seurat::DotPlot(mbrain3, features = paste0("precursor_", 2*(i-1)+c(2,1)),
                           group.by = "celltype") + ggplot2::theme_classic()
  plot3 <- plot3 + ggplot2::scale_x_discrete(labels = c("RNA expression", "ATAC expression (Sum)"))
  
  p <- cowplot::plot_grid(plot1|plot2, plot3, ncol = 1, nrow = 2)
  cowplot::save_plot(filename = paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_precursor_", gene, ".png"), 
                     p, base_asp = 1.2, base_height = 8, device = "png")
}

##################################

mat_x <- as.matrix(Matrix::t(mbrain3[["ATAC"]]@data))
mat_y <- as.matrix(Matrix::t(mbrain3[["SCT"]]@data))

# genes <- rownames(de_combined[[1]][1:50,])
genes <- unique(chiyun$gene)
genes <- genes[which(genes %in% colnames(mat_y))]

precursor_score_normalized <- lapply(genes, function(gene){
  idx <- which(mbrain3[["ATAC"]]@links$gene == gene)
  peak_names <- mbrain3[["ATAC"]]@links$peak[idx]
  tmp <- mat_x[,which(colnames(mat_x) %in% peak_names), drop = F]
  tmp <- sapply(1:ncol(tmp), function(j){
    tmp[,j]/max(tmp[,j])
  })
  vec1 <- rowSums(tmp)
  vec2 <- mat_y[,which(colnames(mat_y) == gene)]
  
  vec1 <- pmin(vec1, stats::quantile(vec1, probs = 0.9))
  # vec2 <- pmin(vec2, stats::quantile(vec2, probs = 0.9))
  if(max(vec1) > 0) vec1 <- vec1/max(vec1)
  if(max(vec2) > 0) vec2 <- vec2/max(vec2)
  
  tmp <- cbind(vec1, vec2)
  colnames(tmp) <- c("precursor", "rna")
  tmp
})

total_peak <- rowSums(sapply(precursor_score_normalized, function(x){
  x[,1]
}))
total_rna <- rowSums(sapply(precursor_score_normalized, function(x){
  x[,2]
}))

precursor_mat <- cbind(total_peak, total_rna, log(total_peak))
colnames(precursor_mat) <- paste0("precursor", 1:ncol(precursor_mat))
mbrain3[["precursor"]] <- Seurat::CreateDimReducObject(embedding = precursor_mat, 
                                                       key = "precursor_", assay = "RNA")

plot1 <- Seurat::FeaturePlot(mbrain3, features = "precursor_2", reduction = "wnn.umap")
plot1 <- plot1 + ggplot2::labs(title = "DE RNA", x = "wnnUMAP_1", y = "wnnUMAP_2")

plot2 <- Seurat::FeaturePlot(mbrain3, features = "precursor_1", reduction = "wnn.umap")
plot2 <- plot2 + ggplot2::labs(title = "Linked ATAC (Precursor score)", x = "wnnUMAP_1", y = "wnnUMAP_2")

plot3 <- Seurat::DotPlot(mbrain3, features = c("precursor_2", "precursor_3"),
                         group.by = "celltype") + ggplot2::theme_classic()
plot3 <- plot3 + ggplot2::scale_x_discrete(labels = c("RNA expression", "Precursor score"))

p <- cowplot::plot_grid(plot1|plot2, plot3, ncol = 1, nrow = 2)
cowplot::save_plot(filename = paste0("../../../../out/figures/Writeup3c/Writeup3c_10x_embryo_precursor_all_oligo.png"), 
                   p, base_asp = 1.2, base_height = 8, device = "png")
