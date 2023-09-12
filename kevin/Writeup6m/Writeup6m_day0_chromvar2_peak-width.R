rm(list=ls())
library(Seurat); library(Signac)
library(GenomicRanges); library(GenomeInfoDb); library(IRanges)
library(JASPAR2020); library(TFBSTools); library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract.RData")

Seurat::DefaultAssay(all_data) <- "ATAC"
mf <- all_data[["ATAC"]]@meta.features["sequence.length"]

p1 <- ggplot2::ggplot(mf, aes(x=sequence.length)) + ggplot2::geom_histogram()
p1 <- p1 + ggplot2::xlab("Peak sequence length")

ggplot2::ggsave(filename = "../../../../out/figures/Writeup6m/Writeup6m_peak-width.png",
                p1, device = "png", width = 13, height = 6, units = "in")


######

zz <- all_data[["ATAC"]]
zz2 <- zz@ranges
zz3 <- zz2@ranges
