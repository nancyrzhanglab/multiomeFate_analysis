rm(list=ls())
library(Seurat); library(Signac)
library(GenomicRanges); library(GenomeInfoDb); library(IRanges)
library(JASPAR2020); library(TFBSTools); library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

motif_focus <- "FOS"
load(paste0("../../../../out/kevin/Writeup6m/Writeup6m_day0-vs-COCL2week5_chromvar2_", motif_focus, ".RData"))

n <- nrow(chromvar_mat)
day0_idx <- grep("day0", rownames(chromvar_mat))
week5_idx <- grep("week5", rownames(chromvar_mat))
status_vec <- rep(NA, n)
status_vec[day0_idx] <- "day0"
status_vec[week5_idx] <- "week5"
status_vec <- factor(status_vec)
names(status_vec) <- rownames(chromvar_mat)

###############

pdf(paste0("../../../../out/figures/Writeup6m/Writeup6m_day0-vs-COCL2week5_chromvar2_", motif_focus, "_violin.pdf"), 
    onefile = T, width = 5, height = 5)

for(var_idx in 1:ncol(chromvar_mat)){
  df <- data.frame(status = status_vec,
                   score = chromvar_mat[,var_idx])
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=status, y=score))
  p1 <- p1 + ggplot2::geom_violin(trim=T, scale = "width")
  p1 <- p1 + ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha = 0.3, size = 0.5)
  p1 <- p1 + Seurat::NoLegend()
  p1 <- p1 + ggplot2::xlab("Winner/Loser status")
  p1 <- p1 + ggplot2::geom_boxplot(width=0.05)
  p1 <- p1 + ggplot2::stat_summary(fun=mean, geom="point", shape=16, size=3, color="red")
  p1 <- p1 + ggplot2::ggtitle(colnames(chromvar_mat)[var_idx])
  print(p1)
}

dev.off()

##########################################

# let's try a dimension-reduction
seurat_obj <- Seurat::CreateSeuratObject(counts = t(chromvar_mat), 
                                         meta.data = data.frame(status_vec))
seurat_obj[["RNA"]]@var.features <- colnames(chromvar_mat)
seurat_obj <- Seurat::ScaleData(seurat_obj)
set.seed(10)
seurat_obj <- Seurat::RunPCA(seurat_obj, verbose = F)
seurat_obj <- Seurat::RunUMAP(seurat_obj, verbose = F, dims = 1:5)

plot1 <-Seurat::DimPlot(seurat_obj, reduction = "umap",
                        group.by = "status_vec", label = TRUE,
                        repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Day0 and COCL week5 cells:\nChromvar2 scores for ", motif_focus, " motifs"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup6m/Writeup6m_day0-vs-COCL2week5_chromvar2_", motif_focus, "_umap.png"), 
                plot1, device = "png", width = 7, height = 5, units = "in")





