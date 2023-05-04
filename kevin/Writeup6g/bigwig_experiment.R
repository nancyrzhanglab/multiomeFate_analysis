# the JUN gene
library(rtracklayer); library(GenomicRanges)
which_range <- GenomicRanges::GRanges("chr1", IRanges(58780000, 58785000))
zz <- rtracklayer::import("EGFR_BC.bw", which = which_range)
head(zz)
head(zz@elementMetadata)

x_vec <- zz@ranges@start
y_vec <- zz@elementMetadata$score

png("~/project/Multiome_fate/out/figures/Writeup6g/bigwig_test.png",
    width = 2500, height = 1200, units = "px", res = 300)
plot(NA, xlim = range(x_vec), ylim = range(y_vec), 
     xlab = "Basepair", ylab = "Score")
polygon(x = c(x_vec[1], x_vec, x_vec[length(x_vec)]),
        y = c(0, y_vec, 0),
        col = "gray")
graphics.off()

#########################

zz <- wig::import_wig("EGFR_NGFR_BC.wig")
head(zz)