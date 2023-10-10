library(Seurat)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# ==============================================================================
# Read data
# ==============================================================================

## Read Day10 growth potentials
tp_early <- "day10"
treatment <- "CIS"
load(paste0("/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cis_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "COCL2"
load(paste0("/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cocl2_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "DABTRAM"
load(paste0("/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/Multiome_fate/out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
dabtram_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]

## Read chromVars
load('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/chromVar_day10_data.RData')
cis_day10_chromVar_mat <- t(chromvar_results_cis[,names(cis_d10_imputed)])
cocl2_day10_chromVar_mat <- t(chromvar_results_cocl2[,names(cocl2_d10_imputed)])
dabtram_day10_chromVar_mat <- t(chromvar_results_dabtram[,names(dabtram_d10_imputed)])

## Read Day10 growth potentials
motifs <- read.csv('/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/data/motif_info.csv')

# ==============================================================================
# Calculate correlation data
# ==============================================================================
## CIS
cis_cor_vec <- sapply(1:ncol(cis_day10_chromVar_mat), function(j){
  res <- stats::cor.test(cis_d10_imputed, cis_day10_chromVar_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
cis_cor_vec <- as.data.frame(t(cis_cor_vec))
colnames(cis_cor_vec) <- c("correlation", "p.value")
rownames(cis_cor_vec) <- colnames(cis_day10_chromVar_mat)
cis_cor_vec$motif_code <- rownames(cis_cor_vec)
cis_cor_vec <- merge(cis_cor_vec, motifs, by='motif_code')
rownames(cis_cor_vec) <- cis_cor_vec$motif_names
cis_cor_vec <- cis_cor_vec[, c("correlation", "p.value")]

## COCL2
cocl2_cor_vec <- sapply(1:ncol(cocl2_day10_chromVar_mat), function(j){
  res <- stats::cor.test(cocl2_d10_imputed, cocl2_day10_chromVar_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
cocl2_cor_vec <- as.data.frame(t(cocl2_cor_vec))
colnames(cocl2_cor_vec) <- c("correlation", "p.value")
rownames(cocl2_cor_vec) <- colnames(cocl2_day10_chromVar_mat)
cocl2_cor_vec$motif_code <- rownames(cocl2_cor_vec)
cocl2_cor_vec <- merge(cocl2_cor_vec, motifs, by='motif_code')
rownames(cocl2_cor_vec) <- cocl2_cor_vec$motif_names
cocl2_cor_vec <- cocl2_cor_vec[, c("correlation", "p.value")]

## DABTRAM
dabtram_cor_vec <- sapply(1:ncol(dabtram_day10_chromVar_mat), function(j){
  res <- stats::cor.test(dabtram_d10_imputed, dabtram_day10_chromVar_mat[,j],
                         alternative = "two.sided",
                         method = "pearson")
  c(res$estimate, res$p.value)
})
dabtram_cor_vec <- as.data.frame(t(dabtram_cor_vec))
colnames(dabtram_cor_vec) <- c("correlation", "p.value")
rownames(dabtram_cor_vec) <- colnames(dabtram_day10_chromVar_mat)
dabtram_cor_vec$motif_code <- rownames(dabtram_cor_vec)
dabtram_cor_vec <- merge(dabtram_cor_vec, motifs, by='motif_code')
rownames(dabtram_cor_vec) <- dabtram_cor_vec$motif_names
dabtram_cor_vec <- dabtram_cor_vec[, c("correlation", "p.value")]

# ==============================================================================
# Saving
# ==============================================================================
correlation_list <- c(cis_cor_vec, cocl2_cor_vec, dabtram_cor_vec)
save(date_of_run, session_info,
     correlation_list,
     file = "/home/mnt/weka/nzh/team/emiliac/nzhanglab/project/emiliac/lineage_trace/task_identify_features_predictive_growth_and_lineage_specific/outs/day10_chromVar_day10_growth_potential_for_week5_correlation.RData")