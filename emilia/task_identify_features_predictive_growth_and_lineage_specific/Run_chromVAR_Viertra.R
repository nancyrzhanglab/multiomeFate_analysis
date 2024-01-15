library(Seurat)
library(GenomicRanges)

# ==============================================================================
# Read data
# ==============================================================================
day0 <- readRDS("/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/Raw_and_Processed/day0_with_motifs.rds")
pwms <- readRDS("/home/jingyaq/Minn/resources/Motifs/Vierstra/V2/Vierstra_Archetype_Motifs_V2.0_pwms.rds")