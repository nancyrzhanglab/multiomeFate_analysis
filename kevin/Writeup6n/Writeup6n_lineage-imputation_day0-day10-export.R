rm(list=ls())
library(Seurat)
library(Signac)
library(GenomicRanges)
library(multiomeFate)
library(IRanges)

tp_early <- "day0"
treatment <- "CIS"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cis_d0_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "COCL2"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cocl2_d0_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
dabtram_d0_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]

tp_early <- "day10"
treatment <- "CIS"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cis_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "COCL2"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
cocl2_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]
treatment <- "DABTRAM"
load(paste0("../../../../out/kevin/Writeup6n/Writeup6n_", treatment, "_", tp_early, "_lineage-imputation_stepdown_concise-postprocessed.RData"))
dabtram_d10_imputed <- cell_imputed_count[!is.na(cell_imputed_count)]

load("../../../../out/kevin/Writeup6l/Writeup6l_chromVar_rna-chromvar_lightweight.RData")

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

imputed_list <- list(
  day0_CIS = cis_d0_imputed,
  day0_COCL2 = cocl2_d0_imputed,
  day0_DABTRAM = dabtram_d0_imputed,
  day10_CIS = cis_d10_imputed,
  day10_COCL2 = cocl2_d10_imputed,
  day10_DABTRAM = dabtram_d10_imputed
)

########

correlation_list <- lapply(imputed_list, function(vec){
  rna_mat <- t(all_data[["Saver"]]@data[,names(vec)])
  cor_mat <- sapply(1:ncol(rna_mat), function(j){
    res <- stats::cor.test(vec, rna_mat[,j],
                           alternative = "two.sided",
                           method = "pearson")
    c(res$estimate, res$p.value)
  })
  cor_mat <- t(cor_mat)
  colnames(cor_mat) <- c("correlation", "p.value")
  rownames(cor_mat) <- colnames(rna_mat)
  cor_mat[is.na(cor_mat)] <- 0
  
  cor_mat
})

save(date_of_run, session_info,
     correlation_list,
     file = "../../../../out/kevin/Writeup6n/Writeup6n_lineage-imputation_day0-day10-export.RData")
