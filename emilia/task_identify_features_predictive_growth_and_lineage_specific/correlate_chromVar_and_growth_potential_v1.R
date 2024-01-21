library(Seurat)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================

in_dir <- "/Users/emiliac/Dropbox/Thesis/Lineage_trace/data/Shaffer_lab/"
out_dir <- "/Users/emiliac/Dropbox/Thesis/Lineage_trace/outputs/task4_identify_genes_corr_growth_and_lineage_specific/"
## Read Day10 growth potentials
tp_early <- "day10"
treatment <- "CIS"
load(paste0(in_dir, 'Growth_potential/Writeup6r_', treatment, "_", tp_early, "_lineage-imputation_postprocess.RData"))
cis_d10_imputed <- cell_imputed_score[!is.na(cell_imputed_score)]
treatment <- "COCL2"
load(paste0(in_dir, "Growth_potential/Writeup6r_", treatment, "_", tp_early, "_lineage-imputation_postprocess.RData"))
cocl2_d10_imputed <- cell_imputed_score[!is.na(cell_imputed_score)]
treatment <- "DABTRAM"
load(paste0(in_dir, "Growth_potential/Writeup6r_", treatment, "_", tp_early, "_lineage-imputation_postprocess.RData"))
dabtram_d10_imputed <- cell_imputed_score[!is.na(cell_imputed_score)]

## Read chromVars
cis_day10_chromVar_mat <- readRDS(paste0(in_dir, 'ChromVar/Vierstra/mat.motifs_day10_CIS_Vierstra_annotations.rds'))
cis_day10_chromVar_mat <- t(cis_day10_chromVar_mat@assays@data@listData[["z"]])
cocl2_day10_chromVar_mat <- readRDS(paste0(in_dir, 'ChromVar/Vierstra/mat.motifs_day10_COCL2_Vierstra_annotations.rds'))
cocl2_day10_chromVar_mat <- t(cocl2_day10_chromVar_mat@assays@data@listData[["z"]])
dabtram_day10_chromVar_mat <- readRDS(paste0(in_dir, 'ChromVar/Vierstra/mat.motifs_day10_DABTRAM_Vierstra_annotations.rds'))
dabtram_day10_chromVar_mat <- t(dabtram_day10_chromVar_mat@assays@data@listData[["z"]])


# ==============================================================================
# Calculate correlation data
# ==============================================================================
## CIS
cis_day10_chromVar_mat <- cis_day10_chromVar_mat[names(cis_d10_imputed), ]
print(all(rownames(cis_day10_chromVar_mat) == names(cis_d10_imputed)))
cis_cor_vec <- sapply(1:ncol(cis_day10_chromVar_mat), function(j){
  res <- stats::cor.test(cis_d10_imputed, cis_day10_chromVar_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
cis_cor_vec <- as.data.frame(t(cis_cor_vec))
colnames(cis_cor_vec) <- c("correlation", "p.value")
rownames(cis_cor_vec) <- colnames(cis_day10_chromVar_mat)

## COCL2
cocl2_day10_chromVar_mat <- cocl2_day10_chromVar_mat[names(cocl2_d10_imputed), ]
print(all(rownames(cocl2_day10_chromVar_mat) == names(cocl2_d10_imputed)))
cocl2_cor_vec <- sapply(1:ncol(cocl2_day10_chromVar_mat), function(j){
  res <- stats::cor.test(cocl2_d10_imputed, cocl2_day10_chromVar_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
cocl2_cor_vec <- as.data.frame(t(cocl2_cor_vec))
colnames(cocl2_cor_vec) <- c("correlation", "p.value")
rownames(cocl2_cor_vec) <- colnames(cocl2_day10_chromVar_mat)


## DABTRAM
dabtram_day10_chromVar_mat <- dabtram_day10_chromVar_mat[names(dabtram_d10_imputed), ]
print(all(rownames(dabtram_day10_chromVar_mat) == names(dabtram_d10_imputed)))
dabtram_cor_vec <- sapply(1:ncol(dabtram_day10_chromVar_mat), function(j){
  res <- stats::cor.test(dabtram_d10_imputed, dabtram_day10_chromVar_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
dabtram_cor_vec <- as.data.frame(t(dabtram_cor_vec))
colnames(dabtram_cor_vec) <- c("correlation", "p.value")
rownames(dabtram_cor_vec) <- colnames(dabtram_day10_chromVar_mat)


# ==============================================================================
# Saving
# ==============================================================================
correlation_list <- list(cis_cor_vec, cocl2_cor_vec, dabtram_cor_vec)
names(correlation_list) <- c('cis_cor_vec', 'cocl2_cor_vec', 'dabtram_cor_vec')

save(date_of_run, session_info,
     correlation_list,
     file = paste0(out_dir, "day10_chromVar_day10_growth_potential_for_week5_correlation_writeup6r.RData"))


