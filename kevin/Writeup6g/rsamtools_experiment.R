library("Rsamtools")

bamFile <- Rsamtools::BamFile("EGFR_Merged_Reads.bam")
bamFile
GenomeInfoDb::seqinfo(bamFile)

# read in JUN
gr <- GenomicRanges::GRanges(seqnames = "chr1",
              ranges = IRanges(start = c(58780000), end = c(58785000)))
params <- Rsamtools::ScanBamParam(which = gr, what = Rsamtools::scanBamWhat())
aln <- Rsamtools::scanBam(bamFile, param = params)
class(aln)
names(aln)
length(aln)
names(aln[[1]])
head(aln[[1]]$pos)
quantile(aln[[1]]$qwidth)
quantile(aln[[1]]$mapq)
head(aln[[1]]$flag)

# now for a pileup
pileupParam <- Rsamtools::PileupParam(max_depth=5000,
                                      distinguish_strands=F, 
                                      distinguish_nucleotides=F)
pileup_res <- Rsamtools::pileup(bamFile, scanBamParam=params,
               pileupParam=pileupParam)
class(pileup_res)
head(pileup_res)
dim(pileup_res)

x_vec <- pileup_res$pos
y_vec <- pileup_res$count
png("~/project/Multiome_fate/out/figures/Writeup6g/pileup_test.png",
    width = 2500, height = 1200, units = "px", res = 300)
plot(NA, xlim = range(x_vec), ylim = range(y_vec), 
     xlab = "Basepair", ylab = "Score")
polygon(x = c(x_vec[1], x_vec, x_vec[length(x_vec)]),
        y = c(0, y_vec, 0),
        col = "gray")
graphics.off()


# cleanup
close(bamFile)
yieldSize(bamFile) <- NA

