rm(list=ls())
library(Seurat)
load("../../../../out/kevin/Writeup4a/2022-02-09_time0_formatted.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

plot1 <- Seurat::VlnPlot(t0_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_QC_metric.png"),
                plot1, device = "png", width = 10, height = 6, units = "in")

plot1 <- Seurat::FeatureScatter(t0_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- plot1 + Seurat::NoLegend()
plot2 <- Seurat::FeatureScatter(t0_obj, feature1 = "nCount_Lineage", feature2 = "nFeature_Lineage")
plot2 <- plot2 + Seurat::NoLegend()
plot1 <- plot1 + plot2
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_QC_scatterplot.png"),
                plot1, device = "png", width = 10, height = 6, units = "in")

################

tmp <- t0_obj[["Lineage"]]@counts; tmp@x <- rep(1, length(tmp@x))
colsum_vec <- Matrix::colSums(tmp)
png("../../../../out/figures/Writeup4a/2022-02-09_lineage_histogram.png", 
    height = 900, width = 1500, res = 300, units = "px")
graphics::hist(colsum_vec, breaks = seq(min(colsum_vec)-.5, max(colsum_vec)+.5, by = 1),
               col = "gray", main = "Number of lineage barcodes for a cell",
               xlab = "Number of lineages", ylab = "Frequency")
graphics.off()

##################

set.seed(10)
tmp <- Seurat::SCTransform(t0_obj)

jackpot_genes <- c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "GADPH", "CCNA2", "MMP1",
                   "IER2", "MYC", "EEF2", "SERTAD4", "IFIT1", "PMAIP1", "DDX58", "CYR61", "APCDD1", "SOX3", "TEAD", "TEAD4", "ATF4", "BATF") # others
jackpot_genes <-  intersect(rownames(t0_obj[["RNA"]]), jackpot_genes)
gene_vec <- unique(c(Seurat::VariableFeatures(tmp, assay = "SCT"),
                     jackpot_genes))
set.seed(10)
t0_obj <- Seurat::SCTransform(t0_obj, 
                              residual.features = gene_vec)
t0_obj <- Seurat::RunPCA(t0_obj, npcs = 30, verbose = FALSE)
set.seed(10)
t0_obj <- Seurat::RunUMAP(t0_obj, reduction = "pca", dims = 1:30)
t0_obj <- Seurat::FindNeighbors(t0_obj, reduction = "pca", dims = 1:30, verbose = FALSE)
t0_obj <- Seurat::FindClusters(t0_obj, resolution = 1, verbose = FALSE)

plot1 <- Seurat::DimPlot(t0_obj, 
                         reduction = "umap", 
                         group.by = "seurat_clusters", 
                         label = TRUE,
                         repel = TRUE)
plot2 <- Seurat::FeaturePlot(t0_obj, 
                             reduction = "umap", 
                             features = "nCount_RNA")
lineage_cutoff <- t0_obj$nCount_Lineage
lineage_cutoff[lineage_cutoff > 1] <- 2
t0_obj$Lineage_cutoff <- lineage_cutoff
plot3 <- Seurat::FeaturePlot(t0_obj, 
                             reduction = "umap", 
                             features = "Lineage_cutoff")
plot1 <- plot1 + plot2 + plot3
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_umap.png"),
                plot1, device = "png", width = 15, height = 4.5, units = "in")

plot1 <- Seurat::DimPlot(t0_obj, 
                          reduction = "umap", 
                          group.by = "Phase", 
                          label = TRUE,
                          repel = TRUE)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_umap_phase.png"),
                plot1, device = "png", width = 5, height = 4.5, units = "in")

###############

mito_idx <- grep("^MT-", rownames(t0_obj[["SCT"]]@scale.data))
rownames(t0_obj[["SCT"]]@scale.data)[mito_idx]
mito_l2 <- sqrt(Matrix::colSums(t0_obj[["SCT"]]@scale.data[mito_idx,]^2))
all_l2 <- sqrt(Matrix::colSums(t0_obj[["SCT"]]@scale.data^2))
ratio_vec <- mito_l2/all_l2
t0_obj$mt.ratio.sctransform <- ratio_vec
plot1 <- Seurat::FeaturePlot(t0_obj, 
                             reduction = "umap", 
                             features = c("percent.mt", "mt.ratio.sctransform"),
                             ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_umap_mitochondria.png"),
                plot1, device = "png", width = 10, height = 4.5, units = "in")

mito_idx <- grep("^MT-", rownames(t0_obj[["RNA"]]@counts))
Matrix::rowSums(t0_obj[["RNA"]]@counts[mito_idx,])
mito_sum <- Matrix::colSums(t0_obj[["RNA"]]@counts[mito_idx,])
t0_obj$nCount_RNA_mito <- mito_sum

Seurat::Idents(t0_obj) <- rep("1", nrow(t0_obj[["RNA"]]))
plot1 <- Seurat::FeatureScatter(t0_obj, 
                                feature1 = "nCount_RNA_mito", feature2 = "nCount_RNA",
                                group.by = NULL)
plot1 <- plot1 + ggplot2::theme_gray() + Seurat::NoLegend() + ggplot2::coord_fixed()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_QC_scatterplot_mito.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

############

table(Matrix::colSums(t0_obj[["Lineage"]]@counts))
table(Matrix::rowSums(t0_obj[["Lineage"]]@counts))
res <- barcode_assignment(seurat_obj = t0_obj)
t0_obj <- res$seurat_obj
lineage_preprocessing_outs <- res
lineage_preprocessing_outs$seurat_obj <- NULL

save(t0_obj, lineage_preprocessing_outs, date_of_run, session_info,
     file = "../../../../out/kevin/Writeup4a/2022-02-09_time0_preprocessed.RData")

plot1 <- plot_barcode_threshold(res$threshold_df)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_barcode_threshold.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")
plot1 <- plot_lineage_barcodecounts(t0_obj)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_lineage_count.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

###################

col_vec <- scales::hue_pal()(14)
names(col_vec) <- as.character(0:13)
col_vec <- col_vec[c("13","9","4","8","10","3","1","6","5","2","0","7","12","11")] # THIS IS HARD-CODED
plot1 <- Seurat::VlnPlot(t0_obj, 
                         features = "mt.ratio.sctransform", 
                         group.by = "seurat_clusters", 
                         sort = TRUE,
                         pt.size = 0.1, 
                         cols = col_vec)
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_mitochondria_violin.png"),
                plot1, device = "png", width = 8, height = 5, units = "in")

jackpot_genes1 <- c("SOX10", "NGFR", "FN1", "AXL", "EGFR", "NT5E")
jackpot_genes2 <- c("C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                    "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                    "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                    "RUNX2", "LOXL2", "JUN", "PDGFRC", "GADPH", "CCNA2", "MMP1")
other_genes <- c("IER2", "MYC", "EEF2", "SERTAD4", "IFIT1", "PMAIP1", "DDX58", 
                 "CYR61", "APCDD1", "SOX3", "TEAD", "TEAD4", "ATF4", "BATF")
jackpot_genes1 <-  intersect(Seurat::VariableFeatures(t0_obj, assay = "SCT"), jackpot_genes1)
jackpot_genes2 <-  intersect(Seurat::VariableFeatures(t0_obj, assay = "SCT"), jackpot_genes2)
other_genes <-  intersect(Seurat::VariableFeatures(t0_obj, assay = "SCT"), other_genes)

Seurat::DefaultAssay(t0_obj) <- "SCT"
plot1 <- Seurat::FeaturePlot(t0_obj, 
                             slot = "scale.data",
                             reduction = "umap",
                             features = sort(jackpot_genes1),
                             ncol = 3)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_jackpot1.png"),
                plot1, device = "png", width = 12, height = 7, units = "in")
Seurat::DefaultAssay(t0_obj) <- "SCT"
plot1 <- Seurat::FeaturePlot(t0_obj, 
                             slot = "scale.data",
                             reduction = "umap",
                             features = sort(jackpot_genes2),
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_jackpot2.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")
Seurat::DefaultAssay(t0_obj) <- "SCT"
plot1 <- Seurat::FeaturePlot(t0_obj, 
                             slot = "scale.data",
                             reduction = "umap",
                             features = sort(other_genes),
                             ncol = 4)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_jackpot3.png"),
                plot1, device = "png", width = 18, height = 12, units = "in")

####################

lineages_names <- t0_obj$final_lineage
lineages_names <- lineages_names[!lineages_names %in% c("No barcode", "Still multiple", "Too large Naive")]
lineage_table <- sort(table(lineages_names), decreasing = T)
set.seed(10)
lineage_list <- lapply(2:4, function(lineage_size){
  tmp <- names(lineage_table)[which(lineage_table == lineage_size)]
  sample(tmp, 25)
})
names(lineage_list) <- c(2:4)
for(i in 1:length(lineage_list)){
  lis <- lineage_list[[i]]
  cell_name <- rownames(t0_obj@meta.data)[which(t0_obj$final_lineage == lis[1])]
  plot1 <- Seurat::DimPlot(t0_obj, 
                              reduction = "umap",
                              cells.highlight = cell_name)
  plot1 <- plot1 + Seurat::NoLegend()
  
  for(lineage_name in lis[-1]){
    cell_name <- rownames(t0_obj@meta.data)[which(t0_obj$final_lineage == lineage_name)]
    tmp_plot <- Seurat::DimPlot(t0_obj, 
                                reduction = "umap",
                                cells.highlight = cell_name)
    tmp_plot <- tmp_plot + Seurat::NoLegend()
    plot1 <- plot1 + tmp_plot
  }
  
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_example_lineages_size", names(lineage_list)[i], ".png"),
                  plot1, device = "png", width = 12, height = 10, units = "in")
}

##################

# histogram of the Gini coefficient among all the variable genes (based on the counts)
Seurat::DefaultAssay(t0_obj) <- "SCT"
mat <- t0_obj[["RNA"]]@counts[Seurat::VariableFeatures(t0_obj),]
mat <- mat + 0.5
mat <- mat %*% Matrix::Diagonal(x = 1/Matrix::colSums(mat))
mat <- log1p(as.matrix(mat))
# mat <- scale(mat, center = F, scale = T)
gini_vec <- apply(mat, 1, function(x){
  # formula from https://github.com/AndriSignorell/DescTools/blob/master/R/StatsAndCIs.r#L4700
  n <- length(x)
  x <- sort(x)
  res <- 2 * sum(x * 1:n) / (n*sum(x)) - 1 - (1/n)
  pmax(0,  n/(n-1) * res)
})
set.seed(10)
idx <- rank(gini_vec+stats::runif(length(gini_vec), max = 0.001))

jackpot_genes <- c(jackpot_genes1, jackpot_genes2)
text_vec <- rownames(mat)
gene_idx <- which(rownames(mat) %in% jackpot_genes)
gini_vec2 <- c(gini_vec[-gene_idx], gini_vec[gene_idx])
idx2 <- c(idx[-gene_idx], idx[gene_idx])
labeling_vec <- as.factor(c(rep(0, nrow(mat)-length(gene_idx)), rep(1, length(gene_idx))))
text_vec2 <- c(text_vec[-gene_idx], text_vec[gene_idx])
set.seed(10)
df <- data.frame(Gini_value = gini_vec2, Gene_rank = idx2, 
                 labeling = labeling_vec,
                 text = text_vec2)

plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = Gene_rank, y = Gini_value)) 
plot1 <- plot1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
plot1 <- plot1 + ggplot2::scale_colour_manual(values=c("black", "red"))
plot1 <- plot1 + ggrepel::geom_text_repel(data = subset(df, labeling == 1), 
                                          ggplot2::aes(label = text, color = labeling),
                                          max.overlaps = 20)
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup4a/2022-02-09_gini_jackpot.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

##################

# look at cells with high FN1, and see which lineages they are
# cell_idx <- which(t0_obj[["SCT"]]@scale.data["FN1",] >= 10)
# lineages_names <- t0_obj@meta.data[cell_idx,"final_lineage"]
# lineages_names <- lineages_names[!lineages_names %in% c("No barcode", "Still multiple", "Too large Naive")]
# lineage_idx <- t0_obj@meta.data[,"final_lineage"] %in% lineages_names
# table(table(t0_obj@meta.data[lineage_idx,"final_lineage"]))

cell_names <- unlist(lapply(jackpot_genes, function(j){
  colnames(t0_obj[["SCT"]]@scale.data)[which(t0_obj[["SCT"]]@scale.data[j,] > 10)]
}))
table(table(cell_names)) # let's focus on 2 or more
high_cell_names <- names(table(cell_names))[which(table(cell_names) >= 2)]
table(t0_obj$final_lineage[high_cell_names])
lineages_names <- t0_obj$final_lineage[high_cell_names]
lineages_names <- lineages_names[!lineages_names %in% c("No barcode", "Still multiple", "Too large Naive")]
table(table(t0_obj$final_lineage[t0_obj$final_lineage %in% lineages_names]))
length(which(t0_obj$final_lineage=="Lin15283"))
length(which(t0_obj$final_lineage=="Lin8033"))
round(t0_obj[["SCT"]]@scale.data[jackpot_genes, which(t0_obj$final_lineage=="Lin8033")])
round(t0_obj[["SCT"]]@scale.data[jackpot_genes, which(t0_obj$final_lineage=="Lin15283")])
