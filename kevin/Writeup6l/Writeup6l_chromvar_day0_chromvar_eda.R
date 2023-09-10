rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(GGally)
library(ggplot2)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0_chromvar_lightweight_noATAC.RData")

table(all_data$dataset)
name_vec <- names(data.use) 
motif_idx_vec <- unique(c(grep("JUN", name_vec),
                          grep("FOS", name_vec),
                          grep("NFE2", name_vec),
                          grep("TEAD", name_vec)))

pdf("../../../../out/figures/Writeup6l/Writeup6l_chromvar_day0-v2_ap1motifs-correlation-with-nCount_ATAC.pdf", 
    onefile = T, width = 8, height = 5)

for(motif_idx in motif_idx_vec){
  df <- data.frame(nCount_ATAC = all_data$nCount_ATAC,
                   chromvar = all_data[["chromvar"]]@data[motif_idx,])
  cor_val <- stats::cor(df$nCount_ATAC, df$chromvar)
  
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x=nCount_ATAC, y=chromvar))
  p1 <- p1 + geom_point(alpha = 0.3)
  p1 <- p1 + ggplot2::xlab("nCount_ATAC")
  p1 <- p1 + ggplot2::ylab(paste0("Chromvar score for ", name_vec[motif_idx]))
  p1 <- p1 + ggplot2::ggtitle(paste0("Correlation: ", round(cor_val,2)))
  print(p1)
}

dev.off()
