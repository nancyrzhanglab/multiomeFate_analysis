# https://satijalab.org/seurat/articles/atacseq_integration_vignette.html
# https://carmonalab.github.io/scGate.demo/scGate.ATAC-seq.html
# https://stuartlab.org/signac/1.2.0/articles/pbmc_vignette.html
# https://github.com/stuart-lab/signac/discussions/499?sort=top?sort=top
# https://github.com/stuart-lab/signac/issues/826

rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(IRanges)
library(JASPAR2020)
library(TFBSTools)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

load("../../../../out/kevin/Writeup6l/Writeup6l_day0-atac_extract_lightweight.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)
Seurat::DefaultAssay(all_data) <- "ATAC"

## see https://stuartlab.org/signac/articles/motif_vignette.html
## https://stuartlab.org/signac/articles/data_structures.html#the-motif-class
print("Apply TFBSTools::getMatrixSet")
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# fix the paths
all_data[["ATAC"]]@fragments[[1]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time0/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[2]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_CIS/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[3]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_COCL2/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[4]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_DABTRAM/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[5]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_CIS/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[6]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_COCL2/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[7]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_DABTRAM/outs/atac_fragments.tsv.gz"

head(Signac::Annotation(all_data))
head(Signac::Fragments(all_data)[[1]])

# https://github.com/stuart-lab/signac/issues/486
# print("Apply subset")
# main.chroms <- GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
# keep.peaks <- which(as.character(GenomeInfoDb::seqnames(GenomicRanges::granges(all_data[["ATAC"]]))) %in% main.chroms)
# all_data[["ATAC"]] <- subset(all_data[["ATAC"]], features = rownames(all_data[["ATAC"]])[keep.peaks])

print("Apply Signac::CreateMotifMatrix")
motif.matrix <- Signac::CreateMotifMatrix(
  features = GenomicRanges::granges(all_data[["ATAC"]]),
  pwm = pfm,
  genome = "hg38",
  use.counts = FALSE
)

###########

# let's do Signac::CreateMotifMatrix step by step
# https://github.com/stuart-lab/signac/blob/master/R/preprocessing.R#L128

features = GenomicRanges::granges(all_data[["ATAC"]])
pwm = pfm
genome = "hg38"
score = FALSE
use.counts = FALSE
sep = c("-", "-")

genome <- BSgenome::getBSgenome(genome = genome)
miss_sn <- !(as.character(seqnames(x = features)) %in% seqlevels(x = genome))

# https://github.com/GreenleafLab/motifmatchr/blob/master/R/match_motifs.R
motif_ix <- motifmatchr::matchMotifs(
  pwms = pwm,
  subject = features,
  genome = genome,
  out = "scores"
)