rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

load("../../../../out/kevin/Writeup6j/Writeup6j_DABTRAM_day0_lineage-imputation_stepdown-LOOCV_concise-postprocessed.RData")
dabtram_imputed <- cell_imputed_count
load("../../../../out/kevin/Writeup6j/Writeup6j_COCL2_day0_lineage-imputation_stepdown-LOOCV_concise-postprocessed.RData")
cocl2_imputed <- cell_imputed_count

cell_names <- intersect(names(dabtram_imputed), names(cocl2_imputed))
dabtram_imputed <- dabtram_imputed[cell_names]
cocl2_imputed <- cocl2_imputed[cell_names]
any(is.na(dabtram_imputed)); any(is.na(cocl2_imputed))

# dabtram_imputed isn't log10
dabtram_imputed <- log10(exp(dabtram_imputed))

stats::cor(dabtram_imputed, cocl2_imputed)
stats::cor(10^dabtram_imputed, 10^cocl2_imputed)

df <- data.frame(cocl2_day0 = cocl2_imputed,
                 dabtram_day0 = dabtram_imputed)
p1 <- GGally::ggpairs(df, 
                     lower = list(continuous = GGally::wrap("points", alpha = 0.2, shape = 16)),
                     progress = FALSE) + ggplot2::theme_bw()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup6n/Writeup6n_lineage-imputation_day0-comparison.png",
                p1, device = "png", width = 6, height = 6, units = "in")

#####################

load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_fasttopics_DABTRAM.RData")
topic_dabtram <- topic_res
load("../../../../out/kevin/Writeup4e/Writeup4e_timeAll_fasttopics_COCL2.RData")
topic_cocl2 <- topic_res

