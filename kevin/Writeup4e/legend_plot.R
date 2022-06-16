rm(list=ls())
library(plotfunctions)

pdf("../../out/figures/Writeup4e_pdfs/enrichment_legend.pdf",
    height = 5, width = 5)
plotfunctions::emptyPlot(1,1, main='Test plot\nfor lineage enrichment values', axes=FALSE)
col_vec <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(100)
plotfunctions::gradientLegend(valRange=c(0, 1), 
                              color = col_vec,
                              length = 0.75,
                              n.seg=3, pos=.5, side=1)
graphics.off()