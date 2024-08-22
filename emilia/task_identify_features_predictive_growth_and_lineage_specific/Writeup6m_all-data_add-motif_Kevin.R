rm(list=ls())
library(Seurat); library(Signac)
library(GenomicRanges); library(GenomeInfoDb); library(IRanges)
library(JASPAR2020); library(TFBSTools); library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges); library(IRanges)

load("../../../../out/kevin/Writeup6b/Writeup6b_all-data.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# now to add the chromvar part
Seurat::DefaultAssay(all_data) <- "ATAC"


## see https://stuartlab.org/signac/articles/motif_vignette.html
## https://stuartlab.org/signac/articles/data_structures.html#the-motif-class
print("Apply TFBSTools::getMatrixSet")
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020::JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# # https://github.com/stuart-lab/signac/issues/486
# to resolve the error:
# ````
# Error in .getOneSeqFromBSgenomeMultipleSequences(x, name, start, NA, width,  :
#                                                    sequence GL000194.1 not found
# ````                                                 
print("Apply subset")
main.chroms <- GenomeInfoDb::standardChromosomes(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(GenomeInfoDb::seqnames(GenomicRanges::granges(all_data[["ATAC"]]))) %in% main.chroms)
all_data[["ATAC"]] <- subset(all_data[["ATAC"]], features = rownames(all_data[["ATAC"]])[keep.peaks])

# fix the paths
all_data[["ATAC"]]@fragments[[1]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time0/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[2]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_CIS/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[3]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_COCL2/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[4]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_time10_DABTRAM/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[5]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_CIS/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[6]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_COCL2/outs/atac_fragments.tsv.gz"
all_data[["ATAC"]]@fragments[[7]]@path <- "/home/stat/nzh/team/kevinl1/project/Multiome_fate/BarcodeOutputs/2022_02/Cellranger_count_output/2022_05_19_arc_week5_DABTRAM/outs/atac_fragments.tsv.gz"

## see https://github.com/stuart-lab/signac/blob/master/R/motifs.R
## see https://github.com/stuart-lab/signac/issues/429
print("Apply Signac::CreateMotifMatrix")
motif.matrix <- Signac::CreateMotifMatrix(
  features = GenomicRanges::granges(all_data[["ATAC"]]),
  pwm = pfm,
  genome = "hg38",
  use.counts = FALSE
)

###################################

print("Apply motifmatchr::matchMotifs")
motif.positions <- motifmatchr::matchMotifs(
  pwms = pfm,
  subject = GenomicRanges::granges(all_data[["ATAC"]]),
  out = "positions",
  genome = "hg38"
)

###################################

print("Apply Signac::CreateMotifObject")
motif <- Signac::CreateMotifObject(
  data = motif.matrix,
  positions = motif.positions,
  pwm = pfm
)

print("Apply Seurat::SetAssayData")
all_data[["ATAC"]] <- Seurat::SetAssayData(
  object = all_data[["ATAC"]],
  slot = 'motifs',
  new.data = motif
)

####################################

print("Apply Signac::RegionStats")
Seurat::DefaultAssay(all_data) <- "ATAC"
all_data <- Signac::RegionStats(
  object = all_data, 
  genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
  assay = "ATAC"
)


print("Saving")
save(date_of_run, session_info, 
     all_data, 
     file = "../../../../out/kevin/Writeup6m/Writeup6m_all-data.RData")
