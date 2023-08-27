rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(ggplot2) 

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

tab_mat <- table(all_data$assigned_lineage, all_data$dataset)

# start binning the cells at day0
table(tab_mat[,"day0"])

tab_vec <- tab_mat[,"day0"]
tab_vec <- tab_vec[tab_vec > 0]
bin_vec <- rep(NA, length(tab_vec))
names(bin_vec) <- names(tab_vec)
bin_vec[which(tab_vec == 1)] <- "1"
bin_vec[which(tab_vec == 2)] <- "2"
bin_vec[which(tab_vec == 3)] <- "3"
bin_vec[intersect(which(tab_vec >= 4), which(tab_vec <= 6))] <- "4-6"
bin_vec[intersect(which(tab_vec >= 7), which(tab_vec <= 9))] <- "7-9"
bin_vec[intersect(which(tab_vec >= 10), which(tab_vec <= 12))] <- "10-12"
bin_vec <- factor(bin_vec, levels = c("1", "2", "3", "4-6", "7-9", "10-12"))

cell_names <- colnames(all_data)[which(all_data$dataset == "day0")]
cell_lineages <- all_data$assigned_lineage[cell_names]
cell_bin_vec <- bin_vec[cell_lineages]

###############

name_vec <- names(data.use) 
motif_idx_vec <- unique(c(grep("JUN", name_vec),
                          grep("FOS", name_vec),
                          grep("NFE2", name_vec),
                          grep("TEAD", name_vec)))

################

pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_day0_violin-all-cells_ap1motifs.pdf", 
    onefile = T, width = 8, height = 5)

for(motif_idx in motif_idx_vec){
  df <- data.frame(bin = cell_bin_vec,
                   chromvar = all_data[["chromvar"]]@data[motif_idx,cell_names])
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=cell_bin_vec, y=chromvar))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width")
  p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::xlab("Day0 Lineage size (Binned)")
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  p1 <- p1 + ggplot2::ggtitle(paste0("Chromvar score for ", name_vec[motif_idx]))
  print(p1)
}

dev.off()

#############

# bin by nCount_ATAC and nFeature_ATAC and TSS.enrichment and pct_reads_in_peaks
variable_vec <- c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", 
                  "peak_region_fragments", "pct_reads_in_peaks", 
                  "S.Score", "G2M.Score")

pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_day0_violin-all-cells_metadata-variables.pdf", 
    onefile = T, width = 8, height = 5)

for(variable in variable_vec){
  df <- data.frame(bin = cell_bin_vec,
                   variable = all_data@meta.data[cell_names,variable])
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=cell_bin_vec, y=variable))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width")
  p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::xlab("Day0 Lineage size (Binned)")
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  p1 <- p1 + ggplot2::ggtitle(variable)
  print(p1)
}

dev.off()

#############################


pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_day0_ap1motifs-correlation-with-nCount_ATAC.pdf", 
    onefile = T, width = 8, height = 5)

for(motif_idx in motif_idx_vec){
  df <- data.frame(nCount_ATAC = all_data$nCount_ATAC[cell_names],
                   chromvar = all_data[["chromvar"]]@data[motif_idx,cell_names])
  cor_val <- stats::cor(df$nCount_ATAC, df$chromvar)
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=nCount_ATAC, y=chromvar))
  p1 <- p1 + geom_point(alpha = 0.3)
  p1 <- p1 + ggplot2::xlab("nCount_ATAC")
  p1 <- p1 + ggplot2::ylab(paste0("Chromvar score for ", name_vec[motif_idx]))
  p1 <- p1 + ggplot2::ggtitle(paste0("Correlation: ", round(cor_val,2)))
  print(p1)
}

dev.off()

quantile(Matrix::rowMeans(all_data[["chromvar"]]@data, na.rm = T))
quantile(sparseMatrixStats::rowSds(all_data[["chromvar"]]@data, na.rm = T))
