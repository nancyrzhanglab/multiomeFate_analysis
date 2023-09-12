rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_postprocessed.RData")


day10_forDABTRAM_growth_potential <- all_data$imputed_count
notes <- paste0("Created from Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_postprocess.R",
                " for analyzing the day10_DABTRAM cells based on the lineages sizes at week5_DABTRAM.",
                " Stepdown was used starting with the fasttopic_DABTRAM (for RNA) and the LSI (for ATAC)",
                " components that had a high enough ANOVA score.")

save(day10_forDABTRAM_growth_potential, notes,
     date_of_run, session_info, 
     file = "../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day10_lineage-imputation_stepdown_extract-growth-potential.RData")


