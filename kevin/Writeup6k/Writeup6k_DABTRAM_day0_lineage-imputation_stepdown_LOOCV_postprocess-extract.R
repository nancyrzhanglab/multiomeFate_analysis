# try data fission
rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day0_lineage-imputation_stepdown-LOOCV_postprocessed.RData")

day0_forDABTRAM_growth_potential <- all_data$imputed_count
notes <- paste0("Created from Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_step-down_LOOCV.R",
                " for analyzing the day0 cells based on the lineages sizes at day10_DABTRAM. LOOCV",
                " stepdown was used starting with the fasttopic_DABTRAM (for RNA) and the LSI (for ATAC)",
                " components that had a high enough ANOVA score.")

save(day0_forDABTRAM_growth_potential, notes,
     date_of_run, session_info, 
     file = "../../../../out/kevin/Writeup6k/Writeup6k_DABTRAM_day0_lineage-imputation_stepdown-LOOCV_extract-growth-potential.RData")

