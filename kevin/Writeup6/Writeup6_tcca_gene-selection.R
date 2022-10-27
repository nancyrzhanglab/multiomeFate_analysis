rm(list=ls())

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(tiltedCCA)

load("../../../../out/kevin/Writeup5a/Writeup5a_tcca_RNA-geneActivity.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

####

multiSVD_obj[["common_mat_1"]] <- NULL
multiSVD_obj[["distinct_mat_1"]] <- NULL
multiSVD_obj[["common_dimred_2"]] <- NULL
multiSVD_obj[["distinct_dimred_2"]] <- NULL

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1,
                                                    bool_modality_1_full = T,
                                                    bool_modality_2_full = T)

####

jackpot_genes <- sort(c("SOX10", "MITF", "FN1", "AXL", "EGFR", "NT5E",
                   "C1S", "FRZB", "SERPINB2", "SERPINE1", "NGFR",
                   "SERPINE2", "NDRG1", "FEZF1", "EGR3", "VGF",
                   "WNT5A", "POSTN", "PDGFRB", "NRG1", "VEGFC", "FOSL1",
                   "RUNX2", "LOXL2", "JUN", "PDGFRC", "CD44", "ID3"))

potential_genes <- sort(intersect(colnames(multiSVD_obj$common_mat_1), 
                                  colnames(multiSVD_obj$common_mat_2)))
length(potential_genes)
jackpot_genes <- intersect(jackpot_genes, potential_genes)
rna_mat <- Matrix::t(all_data[["RNA"]]@data[potential_genes,])
geneact_mat <- Matrix::t(all_data[["geneActivity"]]@data[potential_genes,])
n <- nrow(rna_mat)
bool_vec <- sapply(1:length(potential_genes), function(i){
  if(i %% floor(length(potential_genes)/10) == 0) cat('*')
  # print(paste0(i, " out of ", length(potential_genes)))
  bool1 <- length(tiltedCCA:::.nonzero_col(rna_mat, col_idx = i, bool_value = F)) > .1*n
  bool2 <- length(tiltedCCA:::.nonzero_col(geneact_mat, col_idx = i, bool_value = F)) > .1*n
  bool1 & bool2
})
variable_names <- potential_genes[unique(c(which(bool_vec), which(potential_genes %in% jackpot_genes)))]

selection_res <- tiltedCCA:::postprocess_smooth_variable_selection(
  input_obj = multiSVD_obj,
  bool_use_denoised = T,
  bool_include_intercept = T,
  bool_use_metacells = F,
  bool_use_both_modalities = T,
  cor_threshold = 0.95,
  num_variables = 100,
  sd_quantile = 0.05,
  seurat_obj = all_data,
  seurat_assay_1 = "Saver",
  seurat_assay_2 = "geneActivity",
  seurat_slot = "data",
  variable_names = variable_names,
  verbose = 2
)
sort(selection_res$selected_variables)
sort(names(selection_res$alignment_1))

save(date_of_run, session_info, selection_res,
     file = "../../../../out/kevin/Writeup6/Writeup6_tcca_selected-genes.RData")

length(intersect(selection_res$selected_variables, jackpot_genes))
sort(jackpot_genes[jackpot_genes %in% selection_res$selected_variables])
sort(jackpot_genes[!jackpot_genes %in% selection_res$selected_variables])


