rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(GGally)
library(ggplot2)

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")
tab_mat <- table(all_data$assigned_lineage, all_data$dataset)
load("../../../../out/kevin/Writeup6m/Writeup6m_day0_chromvar2_others.RData")

rm_idx <- which(apply(chromvar_mat, 2, function(x){any(is.na(x))}))
if(length(rm_idx) > 0){
  chromvar_mat <- chromvar_mat[,-rm_idx,drop=F]
}

##########

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

pdf("../../../../out/figures/Writeup6m/Writeup6m_chromvar2_day0_violin_day0cells_others.pdf", 
    onefile = T, width = 8, height = 5)

for(j in 1:ncol(chromvar_mat)){
  df <- data.frame(bin = cell_bin_vec,
                   chromvar = chromvar_mat[cell_names,j])
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=cell_bin_vec, y=chromvar))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width")
  p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::xlab("Day0 Lineage size (Binned)")
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  p1 <- p1 + ggplot2::ggtitle(paste0("Chromvar score for ", colnames(chromvar_mat)[j]))
  print(p1)
}

dev.off()

######################################

tab_vec <- tab_mat[,"day10_DABTRAM"]
bin_vec <- rep(NA, length(tab_vec))
names(bin_vec) <- names(tab_vec)
bin_vec[which(tab_vec == 0)] <- "0"
bin_vec[which(tab_vec == 1)] <- "1"
bin_vec[intersect(which(tab_vec >= 2), which(tab_vec <= 5))] <- "2-5"
bin_vec[intersect(which(tab_vec >= 6), which(tab_vec <= 10))] <- "6-10"
bin_vec[intersect(which(tab_vec >= 11), which(tab_vec <= 30))] <- "11-30"
bin_vec[intersect(which(tab_vec >= 31), which(tab_vec <= 150))] <- "31-150"
bin_vec <- factor(bin_vec, levels = c("0", "1", "2-5", "6-10", "11-30", "31-150"))

cell_names <- colnames(all_data)[which(all_data$dataset == "day0")]
cell_lineages <- all_data$assigned_lineage[cell_names]
cell_bin_vec <- bin_vec[cell_lineages]

pdf("../../../../out/figures/Writeup6m/Writeup6m_chromvar2_day0_violin_day0cells_others_byday10DABTRAM.pdf", 
    onefile = T, width = 8, height = 5)

for(j in 1:ncol(chromvar_mat)){
  df <- data.frame(bin = cell_bin_vec,
                   chromvar = chromvar_mat[cell_names,j])
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=cell_bin_vec, y=chromvar))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width")
  p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::xlab("Day0 Lineage size (Binned by Day10 size in DABTRAM)")
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  p1 <- p1 + ggplot2::ggtitle(paste0("Chromvar score for ", colnames(chromvar_mat)[j]))
  print(p1)
}

dev.off()